# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import Molecule
from fragmentation.models import FragSimConf

# from decimal import *
from fragmentation.tests.test_frag_common import FragCommonTests

# from base import tasks


class FragSimConfModelTests(TransactionTestCase):
    def test_file_path(self):
        fs = FragCommonTests.new_frag_sim()
        fsc = fs.frag_sim_conf
        self.assertEqual(
            fsc.file_path("param"),
            FragSimConf.CFM_ID_FOLDER + "param/param_output0.log",
        )
        self.assertEqual(
            fsc.file_path("conf"), FragSimConf.CFM_ID_FOLDER + "conf/param_config.txt"
        )

    def test_file_exist(self):
        fs = FragCommonTests.new_frag_sim()
        fsc = fs.frag_sim_conf
        for fn in FragSimConf.PARAM_FILES:
            self.assertTrue(fsc.file_exist(fn))
            fsc.__setattr__(fn + "_path", "FALSE")
            self.assertFalse(fsc.file_exist(fn))
