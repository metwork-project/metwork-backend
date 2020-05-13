# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import Molecule
from fragmentation.models import FragMolSample
from fragmentation.utils import AdductManager


class AdductTests(TransactionTestCase):
    def test_get_adduct(self):

        PARENT_MASS = 382.32800
        SMILES = "OCC(OC/C=C(C)/CC/C=C(C)/CC/C=C(C)/CC/C=C(C)/C)CO"
        ION_CHARGE = "positive"
        ADDUCT_TARGET = "M+NH4"

        fms = FragMolSample(parent_mass=PARENT_MASS)
        m = Molecule.load_from_smiles(SMILES)
        ad = AdductManager(ION_CHARGE)
        self.assertEqual(ad.get_adduct(m, fms), ADDUCT_TARGET)
