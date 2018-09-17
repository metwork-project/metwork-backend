# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TransactionTestCase
from base.modules import ChemDoodle
import json
from rdkit import Chem

class ChemDoodleTests(TransactionTestCase):

    def test_json_to_mol(self):
        # smiles = 'C/C=C/[C@@](C)(N)c1ccc(O)cc1'
        smiles = 'C/C(F)=C/[C@@](C)(N)c1ccc(O)cc1'
        # smiles = 'C[C@@](N)(/C=C(\F)[NH2])c1ccc([*])cc1'

        json_path = 'base/tests/files/chemdoodle_mol_1.json'
        cd  = ChemDoodle()
        with open(json_path, 'r') as fjson:
            json_str = '[{0}]'.format(fjson.readline())
            json_mol = json.loads(json_str)[0]['m'][0]

        self.assertEqual(cd.json_to_mol(json_mol).smiles(),smiles)

        smarts = '[N,O:1]/[#6](-[#9])=[#6]/[#6@@](-[#6])(-[#7])-[#6]1:[#6]:[#6]:[#6](-[*:2]):[#6]:[#6]:1'
        json_path = 'base/tests/files/chemdoodle_mol_2.json'
        cd  = ChemDoodle()
        with open(json_path, 'r') as fjson:
            json_str = '[{0}]'.format(fjson.readline())
            json_data = json.loads(json_str)
            json_mol = json_data[0]['m'][0]
            json_map = json_data[0]['s']
        mol = cd.json_to_rdkit(json_mol, map=json_map )
        self.assertEqual(Chem.MolToSmarts(mol),smarts)
