# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TransactionTestCase
from base.models import Molecule
from metabolization.modules import ReactionTestManagement
from metabolization.models import Reaction, ReactProcess

class ReactProcessTests(ReactionTestManagement):

    def test_process_reaction_single_reactant(self):
        self.create_reacts(['methylation'])
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

    def test_methods(self):
        self.create_reacts(['methylation'], email='with@file.com')
        self.create_reacts(['methylation_NO_FILE'], email='without@file.com')
        r = Reaction.objects.get(name='methylation')
        r_NO_FILE = Reaction.objects.get(name='methylation_NO_FILE')
        methods = [('reactor',r), ('rdkit',r_NO_FILE)]
        sm = 'OCC'
        expected_smiles = {
            'reactor' : 'COCC',
            'rdkit' : 'CCOCC'}
        reactant = Molecule.load_from_smiles(sm)

        rp = {\
            m[0]: ReactProcess.objects.create(\
                    reaction = m[1], \
                    method = m[0])\
            for m in methods }
        for m in rp:
            rp[m].reactants.add(reactant)
        self.assertFalse(rp['rdkit'].validate())
        r_NO_FILE.smarts = '[N,O,n,o:1]-[H:2]>>[#0:1]-[#6:2]-[#6]'
        r_NO_FILE.save()
        self.assertTrue(rp['rdkit'].validate())
        for m in rp :
            rp[m].run_reaction()
            self.assertTrue(rp[m].achieved)
    # Check if molecule has been created
        expected_mols = \
            { m[0]: Molecule.find_from_smiles(expected_smiles[m[0]]) \
                for m in methods}
        self.assertTrue(not False in expected_mols)
    # Check if reaction_product has been created
        for m in rp:
            self.assertEqual(rp[m].products.count(), 1)
            self.assertTrue(expected_mols[m] in rp[m].products.all())

    def test_process_reaction_double_reactants(self):
        self.create_reacts([ 'diels_alder_cycloaddition'])
        methods = ['reactor', 'rdkit']
        smiles = ['NC1=NC(C=C)=CN1', 'C=CCC']
        smiles = ['C=Cc1c[nH]c(N)n1', 'C=CCC']
        expected_smiles = [\
            'CCC1CCC=C2N=C(N)NC12', \
            'CCC1CC=C2N=C(N)NC2C1']
        ml = [Molecule.load_from_smiles(sm) for sm in smiles]
        r = Reaction.objects.get(name='diels_alder_cycloaddition')
        # r.smarts = '[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4].[#6:6]=[#6:5]>>[#6:1]1[#6:2]=[#6:3][#6:4]-[#6:5]-[#6:6]-1'
        # r.smarts = '[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]1'
        r.smarts = '[#6:1]=,:[#6:2]/[#6:3]=,:[#6:4]/[H].[#6:5]=,:[#6:6]>>[#6:1]1/[#6:2]=,:[#6:3]\[#6:4]-[#6:6]-[#6:5]-1'
        #[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-1
        # r.smarts = '[#6:1]=[#6:2]-[#6:3]=[#6:4].[#6:6]=[#6:5]>>[#6:1]1-[#6:2]=[#6:3]-[#6:4]-[#6:5]-[#6:6]-1'
        r.save()
        for m in methods:
            self.assertTrue(m in r.methods_available())
        for method in methods:
            for i in range(2):
                rp = ReactProcess.objects.create(reaction = r, method = method)
                for m in ml:
                    rp.reactants.add(m)
                rp.save()
                rp.run_reaction()
            # Check if molecule has been created
                expected_mols = {Molecule.find_from_smiles(sm) for sm in expected_smiles}
                self.assertEqual(set(rp.products.all()), expected_mols, \
                    'error for method {0}\n{1}'.format(method, set(rp.products.all())))
                ml.reverse()
