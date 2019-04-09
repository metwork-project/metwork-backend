# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import Molecule
from fragmentation.models import FragSimConf, FragMol
#from decimal import *
from fragmentation.tests.test_frag_common import FragCommonTests
#from base import tasks

class FragSimModelTests(TransactionTestCase):

    def test_frag_mol(self):
        fs = FragCommonTests.new_frag_sim()
        smiles = FragCommonTests.test_data['smiles']
        m = Molecule.load_from_smiles(smiles)
        adduct = 'M+'
        self.assertFalse(fs.get_frag_mol(m, adduct=adduct))
        fs.frag_molecule(m, adduct=adduct)
        self.assertTrue(fs.get_frag_mol(m, adduct=adduct))

