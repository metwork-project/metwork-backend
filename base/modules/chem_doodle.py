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

        m0 = Chem.MolFromSmiles('')
        mw = Chem.RWMol(m0)

        atoms_cd = {}
        bonds_cd = {}


        for atom in json['a']:
            symbol = 'C' if 'l' not in atom else atom['l']
            atoms_cd[atom['i']] = mw.AddAtom(Chem.Atom(symbol))
        for bond in json['b']:
            bond_id = bond['i']
            if 'o' in bond:
                bond_type = self.bond_type(bond['o'])
            else:
                bond_type = self.bond_type(1)
            bonds_cd[bond_id] = bond
            bonds_cd[bond_id]['rd_id'] = mw.AddBond(bond['b'], bond['e'],bond_type) - 1

        atoms_rd = {v: k for k, v in atoms_cd.items()}
        bonds_rd = {v['rd_id']: k for k, v in bonds_cd.items()}

        mol = mw.GetMol()
        Chem.rdmolops.SanitizeMol(mol)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

        # for a in mol.GetAtoms():
        #     a.SetAtomMapNum(a.GetIdx())

        def a_point(atom_rd):
            x, y = ( json['a'][atom_rd][d] for d in ['x','y'] )
            return Point( x,y )

        def mol_atom(idx):
            return mol.GetAtoms()[idx]

    ### Manage stereochemistry

        for bond in mol.GetBonds():
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

    ### Manage chirality
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                s = []
                for b in atom.GetBonds():
                    if 's' in bonds_cd[bonds_rd[b.GetIdx()]]:
                        s.append(bonds_cd[bonds_rd[b.GetIdx()]]['s'])
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

        # print(Chem.MolToSmiles(mol))

        return Molecule.load_from_rdkit(mol)
