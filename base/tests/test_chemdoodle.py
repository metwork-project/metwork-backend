# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TransactionTestCase
from base.modules import ChemDoodle
import json

class ChemDoodleTests(TransactionTestCase):

    def test_json_to_mol(self):
        # smiles = 'C/C=C/[C@@](C)(N)c1ccc(O)cc1'
        smiles = 'C/C(F)=C/[C@@](C)(N)c1ccc(O)cc1'
        json_path = 'base/tests/files/chemdoodle_mol_1.json'
        cd  = ChemDoodle()
        with open(json_path, 'r') as fjson:
            json_str = '[{0}]'.format(fjson.readline())
            json_mol = json.loads(json_str)[0]['m'][0]

        self.assertEqual(cd.json_to_mol(json_mol).smiles(),smiles)
