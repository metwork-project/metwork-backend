# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from rdkit import Chem
import math
# from sympy import intersection
from sympy.geometry import Point, Line, Segment, intersection
from sympy.vector import Vector, CoordSys3D

#Chem.rdDepictor.Compute2DCoords(mr)

class ChemDoodle(object):

    def bond_type(self, value):
        value_dic  = {
            1: Chem.rdchem.BondType.SINGLE,
            2: Chem.rdchem.BondType.DOUBLE,
        }
        return value_dic[value]

    def json_to_mol(self, json):
        from base.models import Molecule
        return Molecule.load_from_rdkit( self.json_to_rdkit(json) )

    def json_to_rdkit(self, json, map={}, mol_type='reactant'):
        from base.models import Molecule

        m0 = Chem.MolFromSmarts('')
        mw = Chem.RWMol(m0)

        atoms_cd = {}
        query_atoms = []

        for atom in json['a']:
            atom_id = atom['i']
            # Query atoms
            if 'q' in atom:
                v = atom['q']['as']['v']
                if 'a' in v:
                    smarts = '[*]'
                    a_defaut = 'C'
                else:
                    smarts = '[{0}]'.format(','.join([a_ for a_ in v ]))
                    a_defaut = v[0]
                query_atoms.append((atom_id,smarts))
                a = Chem.Atom(a_defaut)
            # Non-query atoms
            else:
                symbol = 'C' if 'l' not in atom else atom['l']
                a = Chem.Atom(symbol)
            atoms_cd[atom_id] = mw.AddAtom(a)
        for bond in json['b']:
            bond_id = bond['i']
            if 'o' in bond:
                bond_type = self.bond_type(bond['o'])
            else:
                bond_type = self.bond_type(1)
            mw.AddBond(bond['b'], bond['e'],bond_type) - 1

        bonds_cd1 = {
            i: mw.GetBondBetweenAtoms(b_cd['b'],b_cd['e']).GetIdx()
            for i,b_cd in enumerate(json['b'])
        }

        bonds_cd = {
            '{0}-{1}'.format(b_cd['b'], b_cd['e']): i
            for i,b_cd in enumerate(json['b'])
        }

        atoms_rd = {v: k for k, v in atoms_cd.items()}
        bonds_rd = {v: k for k, v in bonds_cd.items()}

        def a_point(atom_rd):
            x, y = ( json['a'][atom_rd][d] for d in ['x','y'] )
            return Point( x,y )

        def mol_atom(idx):
            return mw.GetAtoms()[idx]

    ### Manage chirality

        for atom in mw.GetAtoms():
            if atom.GetSymbol() == 'C':
                s = []
                for b in atom.GetBonds():
                    id_1 = '{0}-{1}'.format(
                                        b.GetBeginAtomIdx(),
                                        b.GetEndAtomIdx() )
                    if id_1 in bonds_cd:
                        b_cd = json['b'][ bonds_cd[id_1] ]
                    else:
                        b_cd = json['b'][ bonds_cd['{1}-{0}'.format(
                                            b.GetBeginAtomIdx(),
                                            b.GetEndAtomIdx() )] ]
                    if 's' in b_cd:
                        s.append(b_cd['s'])
                    else:
                        s.append('')
                if len(s) == 4 and s.count('protruding')*s.count('recessed') == 1:
                    points = []
                    for i,b in enumerate(atom.GetBonds()):
                        if  b.GetBeginAtomIdx() != atom.GetIdx():
                            a = b.GetBeginAtomIdx()
                        else:
                            a = b.GetEndAtomIdx()

                        if s[i] == 'protruding':
                            z = 1
                        elif s[i] == 'recessed':
                            z = - 1
                        else:
                            z = 0
                        points.append( Point( [i for i in a_point(a)] + [z] ) )
                    int = Line( points[2].midpoint(points[3]), points[1] )\
                            .intersection(
                                Line( points[1].midpoint(points[2]), points[3] ) )[0]
                    N = CoordSys3D('N')
                    coords = [N.i,N.j,N.k]
                    vectors = []
                    for p in points:
                        v = Vector.zero
                        for v_ in [ n*coords[v] for v,n in enumerate(p-int) ]:
                            v += v_
                        vectors.append(v)
                    if vectors[1].cross(vectors[3]).dot(vectors[0]) > 0:
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
                    else:
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)

    ### Manage stereo

        for bond in mw.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                carbons = [ bond.GetBeginAtomIdx(), bond.GetEndAtomIdx() ]
                # Check if stero can be applied
                if True in [ len(mol_atom(c).GetBonds()) > 1 for c in carbons ]:
                    sa = [ \
                            [b.GetBeginAtomIdx() \
                                if b.GetBeginAtomIdx() != c \
                                else b.GetEndAtomIdx()#,
                                for b in mol_atom(c).GetBonds() \
                                    if b.GetBondType() == Chem.rdchem.BondType.SINGLE ] \
                            for c in carbons ]

                    sa = [ sa[i][0] for i in range(2) ]
                    bond.SetStereoAtoms (sa[0], sa[1])

                    if len( intersection(
                            Segment( a_point(sa[0]), a_point(sa[1]) ),
                            Line( a_point(carbons[0]), a_point(carbons[1]) )
                        )) > 0:
                        bond.SetStereo(Chem.rdchem.BondStereo.STEREOTRANS)
                    else:
                        bond.SetStereo(Chem.rdchem.BondStereo.STEREOCIS)

    # replace query_atoms

        if len(query_atoms) > 0:
            for id, smarts in query_atoms:
                m_a = Chem.MolFromSmarts(smarts)
                Chem.SanitizeMol(m_a)
                mw.ReplaceAtom(atoms_cd[id],m_a.GetAtoms()[0])

    # Map atoms
        dic_type = {
            'reactant': 'a1',
            'product': 'a2' }
        i = 1
        for m in map:
            if m['t'] == 'AtomMapping':
                a_rd_id = atoms_cd[ m[ dic_type[mol_type] ] ]
                mw.GetAtoms()[a_rd_id].SetAtomMapNum(i)
                i += 1


        mol = mw.GetMol()
        Chem.rdmolops.SanitizeMol(mol)

        # To keep stereo in MolToSmiles
        Chem.rdDepictor.Compute2DCoords(mol)

        return mol
