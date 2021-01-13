# -*- coding: utf-8 -*-
from __future__ import unicode_literals


from shutil import copyfile
import string
import os
import collections
from metabolization.modules import ReactionTestManagement
from metabolization.models import Reaction, ReactProcess
from base.models import Molecule


class ReactionModelTests(ReactionTestManagement):
    def test_default_name(self):
        default_name = "new_reaction"
        r = Reaction()
        self.assertEqual(r.name, default_name)

    def test_name_unique(self):
        u = self.get_user()
        n = "name_duplicate"
        r_1 = Reaction(name=n, user=u)
        r_1.save()
        r_2 = Reaction(name=n, user=u)
        self.assertRaises(Exception, r_2.save)

    def test_reaction_ready(self):
        reacts = [("methylation", "[N,O:1]>>[*:1]-[#6]"), ("diels_alder", None)]
        ready_on_create = ["methylation"]
        rd = self.create_reacts(reacts)
        for r_name in rd:
            r = rd[r_name]
            self.assertEqual(r.ready(), r_name in ready_on_create)
        r = Reaction.objects.get(name="diels_alder")
        r.smarts = "[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4]-[H].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-1"
        r.save()
        self.assertTrue(r.ready())

    def test_reactants_number(self):
        reacts = [
            ("methylation", "[N,O:1]>>[*:1]-[#6]"),
            (
                "diels_alder",
                "[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4]-[H].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-1",
            ),
            ("error", None),
        ]

        reactants = {
            "methylation": 1,
            "diels_alder": 2,
            "error": 0,
        }
        rd = self.create_reacts(reacts)
        for r_name in rd:
            r = rd[r_name]
            self.assertEqual(r.get_reactants_number(), reactants[r_name])
            self.assertEqual(r.reactants_number, reactants[r_name])

    def test_mass_delta(self):
        r = self.create_reacts([("methylation", "[N,O:1]>>[*:1]-[#6]")])["methylation"]
        m = Molecule.load_from_smiles("CCO")
        self.assertIsNone(r.mass_delta())
        r.run_reaction([m])
        self.assertEqual(round(r.mass_delta(), 6), 14.01565)

    def test_reset_react_process(self):
        # reset react process if not ACTIVE
        r = self.create_reacts(
            [("methylation", "[N,O:1]>>[*:1]-[#6]-[#6]")], activate=False
        )["methylation"]
        co = Molecule.load_from_smiles("CO")
        cco = Molecule.load_from_smiles("CCO")
        self.assertEqual(r.reactprocess_set.count(), 0)
        r.run_reaction([cco])
        self.assertEqual(r.reactprocess_set.count(), 1)
        r.smarts = "[N,O:1]>>[*:1]-[#6]"
        r.save()
        r.run_reaction([co])
        self.assertEqual(r.reactprocess_set.count(), 1)
        r.status_code = Reaction.status.ACTIVE
        r.save()
        r.run_reaction([cco])
        self.assertEqual(r.reactprocess_set.count(), 2)
