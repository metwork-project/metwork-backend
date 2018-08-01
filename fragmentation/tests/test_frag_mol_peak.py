# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import Molecule
from decimal import *
from fragmentation.tests.test_frag_common import FragCommonTests

import hashlib

class FragMilPeakModelTests(TransactionTestCase):

    def test_frag_moleclule(self):
        fs = FragCommonTests.new_frag_sim()
        m = Molecule.load_from_smiles(FragCommonTests.test_data['smiles'])
        fm = fs.frag_molecule(m)
        _expect_res = FragCommonTests.test_data['result']
        _expect_res = [ set([ ( en, val[0], val[1] ) for val in _expect_res[en] ]) for en in _expect_res ]
        expect_res = set([])
        for r in _expect_res:
            expect_res = expect_res | r
        res = set([])
        for fp in fm.fragmolpeak_set.all():
            res = res | {( int(fp.energy), float(fp.mz), float(fp.intensity) )}
        self.assertEqual(res, expect_res)

