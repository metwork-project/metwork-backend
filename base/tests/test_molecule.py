# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from django_rdkit import models
from base.modules.rdkit_functions import RDKit
from base.models import Molecule

class MoleculeModelTests(TransactionTestCase):

    def test_create_mol_from_smiles(self):
        from rdkit import Chem
    # Check create only for valid smiles
        self.assertIsNot(Molecule.create_from_smiles('C(C)(C)(C)(C)'), False)
    # Check not unicity of mol_rdkit field
        self.assertFalse(Molecule.create_from_smiles('C(C)(C)(C)(C)'))
    # Check create only for valid smiles
        self.assertFalse(Molecule.create_from_smiles('C(C)(C)(C)(C)(C)'))
        m = Molecule.create_from_smiles('NC1=NC(C=C)=CN1')
        kekulize_options = {\
            False: 'C=Cc1c[nH]c(N)n1', \
            True: 'C=CC1=C[NH]C(N)=N1'} 
        for k in kekulize_options:
            self.assertEqual(RDKit.mol_to_smiles(m.mol_rdkit, kekulize = k), kekulize_options[k])

    def test_smiles(self):
        sm_init = 'C(C)(C)(C)(C)'
        sm_target = RDKit.mol_to_smiles(RDKit.mol_from_smiles(sm_init))
        m = Molecule.create_from_smiles(sm_init)
        self.assertEqual(m.smiles(), sm_target)

    def test_find_from_smiles(self):
        m = Molecule.create_from_smiles('CC')
        self.assertIsNot(Molecule.find_from_smiles('CC'), False)
        self.assertEqual(Molecule.find_from_smiles('CC'), m)
        self.assertFalse(Molecule.find_from_smiles('C'))

    def test_load_from_smiles(self):
        m = Molecule.create_from_smiles('CCC')
        self.assertIsNot(Molecule.load_from_smiles('CCC'), False)
        self.assertIsNot(Molecule.load_from_smiles('CCCC'), False)
        self.assertFalse(Molecule.load_from_smiles('HHOHH'))

    def load_from_rdkit(self):
        sm = 'CCCC'
        mol_rdkit = RDKit.mol_from_smiles(sm)
        self.assertFalse(Molecule.find_from_smiles(sm))
        m = Molecule.load_from_rdkit(mol_rdkit)
        self.assertTrue(Molecule.find_from_smiles(sm))
        self.assertEqual(m.mol_rdkit, mol_rdkit)
        self.assertEqual(m, Molecule.load_from_rdkit(mol_rdkit))

    def test_major_tautomerization(self):
        smiles = 'NC1=NC2CCCC=C2N1'
        smiles_out = [
            'NC1=NC2=C(CCCC2)N1', 
            'NC1=NC2CCCC=C2N1',
            'NC1=NC2CCCCC2=N1' ]
        mol = Molecule.load_from_smiles(smiles)
        expected_mols = {Molecule.load_from_smiles(sm) for sm in smiles_out}
        self.assertEqual({m for m in mol.major_tautomers()}, expected_mols)

