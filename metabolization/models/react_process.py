# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import subprocess
import time
from django.db import models
from base.models import Molecule
from metabolization.models import Reaction
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from itertools import chain

class ReactProcess(models.Model):

    reaction = models.ForeignKey(Reaction, on_delete=models.CASCADE, default=None, null=True)
    method = models.CharField(max_length=32, default='rdkit')
    method_hash = models.CharField(max_length=32, default='')
    reactants = models.ManyToManyField(Molecule, related_name="reactants", default=None)
    products = models.ManyToManyField(Molecule, related_name="products", default=None)
    achieved = models.BooleanField(default=False)
    status_code = models.SmallIntegerField(\
                        default = 0)

    class status:
        INIT = 0
        READY = 1
        RUNNING = 2
        DONE = 3
        ERROR = 99

    def validate(self):
        methods_available = self.reaction.methods_available()
        if len(methods_available) == 0:
            return False
        if self.method not in methods_available:
            self.method = methods_available[0]
        return \
            ((self.reactants.count() == self.reaction.reactants_number)\
                or (self.reactants.count() == 1))\
            and (self.reactants.count() > 0)\
            and (self.method in methods_available)

    def run_reaction(self):
        if self.validate():
            products = self.run_reaction_(self.method)
            for m in products:
                self.products.add(m)
                self.achieved = True
                self.save()
        return self

    def run_reaction_(self, method):
        reaction = self.reaction
        if method == 'reactor':
            if (self.reactants.count(), self.reaction.reactants_number) == (1,2):
                sm = [self.reactants.all()[0].smiles() for i in range(2)]
            else:
                sm = [m.smiles() for m in self.reactants.all()]
            r_path = reaction.mrv_path()
            products_smiles = []
            for i in range(reaction.reactants_number):
                subproc_list = ['react'] + sm + ['-r', r_path, '--ignore-error']
                products_smiles = products_smiles + \
                    subprocess.check_output(subproc_list)\
                    .decode().split('\n')
                sm.reverse()
            res = {Molecule.load_from_smiles(sm) \
                for sm in products_smiles if sm !=''} - {False}
        elif method == 'rdkit':
            r_smarts = reaction.smarts
            r_smarts = r_smarts.replace('\\','-').replace('/','-')
            rx = rdChemReactions.ReactionFromSmarts(r_smarts)
            products_rdkit = []
            if self.reactants.count() == 1:
                reactant = self.format_reactant(self.reactants.all()[0].mol_rdkit)
                products_rdkit = list(rx.RunReactant(reactant, 0))
            elif self.reactants.count() == 2:
                reactants = [self.format_reactant(m.mol_rdkit) for m in self.reactants.all()]
                for i in range(2):
                    products_rdkit = products_rdkit + list(rx.RunReactants((reactants[i], reactants[1-i])))
            res = {Molecule.load_from_rdkit(m) for m in list(chain(*products_rdkit))} - {False}
        else:
            res = set([])
        self.status_code = ReactProcess.status.DONE
        self.save()
        return res

    def format_reactant(self,mol):
        smarts = self.reaction.smarts
        mol=Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        if 'H' in smarts or '#1' in smarts:
            return Chem.AddHs(mol)
        else:
            return mol

    def wait_run_end(self, timeout = 30):
        begin = time.time()
        while self.status_code != ReactProcess.status.DONE:
            time.sleep(0.5)
            if (time.time() - begin) > timeout:
                print ('\n#### reaction wait_run_end close due to timeout #####\n')
                return self
            else:
                self.refresh_from_db()
        return self
