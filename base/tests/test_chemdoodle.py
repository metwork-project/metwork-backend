# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TransactionTestCase
from base.models import Molecule
from base.modules import ChemDoodle, ChemDoodleJSONError, TestManagement
from metabolization.models import Reaction
import json
from rdkit import Chem
from django.contrib.auth import get_user_model


class ChemDoodleTests(TestManagement):
    def test_json_to_mol(self):
        smiles = "C/C(F)=C/[C@@](C)(N)c1ccc(O)cc1"
        json_path = "base/tests/files/chemdoodle_mol_1.json"
        cd = ChemDoodle()
        with open(json_path, "r") as fjson:
            json_str = "[{0}]".format(fjson.readline())
            json_mol = json.loads(json_str)[0]["m"][0]

        self.assertEqual(cd.json_to_mol(json_mol).smiles(), smiles)

    def test_json_to_smarts(self):
        smarts = "[N,O:1]/[#6](-[#9])=[#6]/[#6@@](-[#6])(-[#7])-[#6]1:[#6]:[#6]:[#6](-[*:2]):[#6]:[#6]:1"
        json_path = "base/tests/files/chemdoodle_mol_2.json"
        cd = ChemDoodle()
        with open(json_path, "r") as fjson:
            json_str = "[{0}]".format(fjson.readline())
            json_data = json.loads(json_str)
            json_mol = json_data[0]["m"][0]
            json_map = json_data[0]["s"]
        mol = cd.json_to_rdkit(json_mol, map=json_map)
        self.assertEqual(Chem.MolToSmarts(mol), smarts)

    def test_json_to_react(self):
        cd = ChemDoodle()

        error_path = {
            "mol_1": "No arrow",
            "reaction_2_lines": "More than one arrow",
            "wrong_direction": "Line in wrong direction",
            "no_product": "No product",
            "too_many_products": "Too many products",
            "ambiguous_mol": "Ambiguous molecule position",
        }
        for file_name, message in error_path.items():
            with open(
                "base/tests/files/chemdoodle_{0}.json".format(file_name), "r"
            ) as fjson:
                json_str = "[{0}]".format(fjson.readline())
                json_data = json.loads(json_str)[0]
            with self.assertRaises(ChemDoodleJSONError) as cm:
                react = cd.json_to_react(json_data)
            self.assertEqual(cm.exception.args[0], message)

        smarts = "[#6:1]-[#6:2](-[N,O,S])=,:[#8:3].[#7:4]-[#6:5]>>[#6:1]-[#6:2](-[#7:4]-[#6:5])=,:[#8:3]"
        json_path = "base/tests/files/chemdoodle_amide_formation_2.json"
        cd = ChemDoodle()

        with open(json_path, "r") as fjson:
            json_str = "[{0}]".format(fjson.readline())
            json_data = json.loads(json_str)[0]
        react = cd.json_to_react(json_data)

        self.assertEqual(react, smarts)
        u = get_user_model().objects.create(email="user@test.com")
        r = Reaction.create_from_smarts(
            smarts=react, user=u, name="test ChemDoodle import"
        )

    def test_mol_to_json(self):
        smiles = "C/C(F)=C/[C@@](C)(N)c1ccc(O)cc1"
        mol = Molecule.load_from_smiles(smiles)
        json_path = "base/tests/files/chemdoodle_mol_1.json"
        cd = ChemDoodle()
        with open(json_path, "r") as fjson:
            json_str = "[{0}]".format(fjson.readline())
            json_mol = json.loads(json_str)[0]["m"][0]
        json_res = cd.mol_to_json(mol)
        self.assertEqual(cd.json_to_mol(json_res).smiles(), smiles)

    def test_react_to_json(self):
        smarts = "[#6:1]-[#6:2](-[N,O,S])=,:[#8:3].[#7:4]-[#6:5]>>[#6:1]-[#6:2](-[#7:4]-[#6:5])=,:[#8:3]"
        json_path = "base/tests/files/chemdoodle_amide_formation_2.json"
        cd = ChemDoodle()
        u = get_user_model().objects.create(email="user@test.com")
        r = Reaction.create_from_smarts(
            smarts=smarts, user=u, name="test ChemDoodle import"
        )

        json_res = cd.react_to_json(r.react_rdkit())
        # print(json.dumps(json_res))
        self.assertEqual(cd.json_to_react(json_res), smarts)
