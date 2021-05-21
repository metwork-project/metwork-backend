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
from base.modules.rdkit_functions import RDKit
from base.models import BaseModel


class ReactProcess(BaseModel):

    reaction = models.ForeignKey(
        Reaction, on_delete=models.CASCADE, default=None, null=True
    )
    reactants = models.ManyToManyField(Molecule, related_name="reactants", default=None)
    products = models.ManyToManyField(Molecule, related_name="products", default=None)
    status_code = models.SmallIntegerField(default=0)

    class status:
        INIT = 0
        READY = 1
        RUNNING = 2
        DONE = 3
        ERROR = 99

    def ready(self):
        return (
            (
                (self.reactants.count() == self.reaction.reactants_number)
                or (self.reactants.count() == 1)
            )
            and (self.reactants.count() > 0)
            and (self.reaction.ready())
        )

    def run_reaction(self):
        if self.ready():
            products = self.run_reaction_()
            for m in products:
                self.products.add(m)
                self.save()
        return self

    def run_reaction_(self):
        reaction = self.reaction
        r_smarts = reaction.smarts
        r_smarts = r_smarts.replace("\\", "-").replace("/", "-")
        rx = rdChemReactions.ReactionFromSmarts(r_smarts)
        products_rdkit = []
        if self.reactants.count() == 1:
            reactant = self.format_reactant(self.reactants.all()[0].mol_rdkit)
            products_rdkit = list(rx.RunReactant(reactant, 0))
        elif self.reactants.count() == 2:
            reactants = [
                self.format_reactant(m.mol_rdkit) for m in self.reactants.all()
            ]
            for i in range(2):
                products_rdkit = products_rdkit + list(
                    rx.RunReactants((reactants[i], reactants[1 - i]))
                )
        res = {Molecule.load_from_rdkit(m) for m in list(chain(*products_rdkit))} - {
            False
        }
        self.status_code = ReactProcess.status.DONE
        self.save()
        return res

    def format_reactant(self, mol):
        smarts = self.reaction.smarts
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        if "H" in smarts or "#1" in smarts:
            mol = Chem.AddHs(mol)
        RDKit.apply_aromaticity(mol)
        return mol

    def wait_run_end(self, timeout=30):
        begin = time.time()
        while self.status_code != ReactProcess.status.DONE:
            time.sleep(0.5)
            if (time.time() - begin) > timeout:
                print("\n#### reaction wait_run_end close due to timeout #####\n")
                return self
            else:
                self.refresh_from_db()
        return self

    def mass_delta(self):
        if self.products.count() > 0:
            rm = sum((r.mass_exact() for r in self.reactants.all()))
            pm = self.products.first().mass_exact()
            return pm - rm
        else:
            return None
