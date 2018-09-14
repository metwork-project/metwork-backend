# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from rdkit import Chem
import math
from sympy import *
from sympy.geometry import *

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
        # double_bonds = []

        for atom in json['a']:
            symbol = 'C' if 'l' not in atom else atom['l']
            atoms_cd[atom['i']] = mw.AddAtom(Chem.Atom(symbol))
        for bond in json['b']:
            bond_id = bond['i']
            if 'o' in bond:
                # if bond['o'] == 2:
                #     double_bonds.append(bond_id)
                bond_type = self.bond_type(bond['o'])
            else:
                bond_type = self.bond_type(1)
            bonds_cd[bond_id] = mw.AddBond(bond['b'], bond['e'],bond_type)

        atoms_rd = {v: k for k, v in atoms_cd.items()}
        bonds_rd = {v: k for k, v in bonds_cd.items()}

        mol = mw.GetMol()
        Chem.rdmolops.SanitizeMol(mol)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    ### Manage stereochemistry
        def a_point(atom_rd):
            x, y = ( json['a'][atom_rd][d] for d in ['x','y'] )
            return Point( x,y )

        def mol_atom(idx):
            return mol.GetAtoms()[idx]

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

                    # print(bond.GetStereo())
                    # print([bond.GetStereoAtoms()[i] for i in range(2)])
                    # print(Chem.MolToSmiles(mol, isomericSmiles=True))
                    # print(mol.GetBonds()[bond.GetIdx()].GetStereo())



        # mol.GetBonds()[2].SetStereoAtoms(0,4)
        # mol.GetBonds()[2].SetStereo(Chem.rdchem.BondStereo.STEREOTRANS)

        # Chem.rdmolops.SanitizeMol(mol)
        print(Chem.MolToSmiles(mol)) #, isomericSmiles=True))

        return Molecule.load_from_rdkit(mol)
