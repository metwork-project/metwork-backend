# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import Molecule
from fragmentation.models import FragSimConf, FragMol
#from decimal import *
from fragmentation.tests.test_frag_common import FragCommonTests
#from base import tasks
from fragmentation.models import FragSimConf
from fragmentation.modules import FragSim


class FragSimModelTests(TransactionTestCase):

    def test_frag_mol(self):
        fs = FragCommonTests.new_frag_sim()
        smiles = FragCommonTests.test_data['smiles']
        m = Molecule.load_from_smiles(smiles)
        adduct = 'M+'
        self.assertFalse(fs.get_frag_mol(m, adduct=adduct))
        fs.frag_molecule(m, adduct=adduct)
        self.assertTrue(fs.get_frag_mol(m, adduct=adduct))

    def test_frag_mol_neg(self):
        fsc = FragSimConf.objects.create(
            param_path = 'param/param_output0_neg.log',
            conf_path =  'conf/param_config_neg.txt',
        )
        fs = FragSim(fsc)
        smiles = "CCC(C)O"
        m = Molecule.load_from_smiles(smiles)
        adduct = 'M-'
        assert not fs.get_frag_mol(m, adduct=adduct)
        fs.frag_molecule(m, adduct=adduct, ion_charge="negative")
        fms = fs.get_frag_mol(m, adduct=adduct)
        assert fms
        spectrum = fms.fragmolspectrum_set.filter(energy=0).first().spectrum
        assert spectrum[0][0] == 17.00328823
