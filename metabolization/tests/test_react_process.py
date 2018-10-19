# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TransactionTestCase
from base.models import Molecule
from metabolization.modules import ReactionTestManagement
from metabolization.models import Reaction, ReactProcess

class ReactProcessTests(ReactionTestManagement):

    def test_process_reaction_single_reactant(self):
        r = self.create_reacts([
            ('methylation', '[N,O:1]>>[*:1]-[#6]') ]) ['methylation']
        sm = 'NC=C(N)C'
        expected_smiles = ['CNC=C(C)N','CNC(C)=CN']
        m = Molecule.load_from_smiles(sm)
        m.save()
        #Reaction.reactions_update()
        r = Reaction.objects.get(name='methylation')
        rp = ReactProcess.objects.create(reaction = r)
        rp.reactants.add(m)
        rp.run_reaction()
    # Check if molecule has been created
        expected_mols = {Molecule.find_from_smiles(sm) for sm in expected_smiles}
        self.assertTrue(not False in expected_mols)
    # Check if reaction_product has been created
        self.assertEqual(rp.products.count(),2)
        for m_exp in expected_mols:
            self.assertTrue(m_exp in rp.products.all())

    def test_process_reaction_double_reactants(self):
        r = self.create_reacts([
            ('diels_alder', None) ]) ['diels_alder']
        smiles = ['C=Cc1c[nH]c(N)n1', 'C=CCC']
        expected_smiles = [\
            'CCC1CCC=C2N=C(N)NC12', \
            'CCC1CC=C2N=C(N)NC2C1']
        ml = [Molecule.load_from_smiles(sm) for sm in smiles]
        smarts = '[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-1'
        r.load_smarts(smarts)
        for i in range(2):
            rp = ReactProcess.objects.create(reaction = r)
            for m in ml:
                rp.reactants.add(m)
            rp.save()
            rp.run_reaction()
        # Check if molecule has been created
            expected_mols = {Molecule.find_from_smiles(sm) for sm in expected_smiles}
            self.assertEqual(set(rp.products.all()), expected_mols)
            ml.reverse()

    def test_mass_delta(self):
        r = self.create_reacts([('methylation', '[N,O:1]>>[*:1]-[#6]')])['methylation']
        m = Molecule.load_from_smiles('CCO')
        rp = r.run_reaction([m])
        self.assertEqual(round(rp.mass_delta(),6), 14.01565)
