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
        m = Molecule.create_from_smiles('CC1=CC=CC=C1')
        kekulize_options = {\
            False: 'Cc1ccccc1', \
            True: 'CC1=CC=CC=C1'}
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
