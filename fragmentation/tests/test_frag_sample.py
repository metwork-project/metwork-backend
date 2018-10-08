# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from django.contrib.auth.models import User
from fragmentation.models import FragSample
from django.contrib.auth import get_user_model
from django.db import IntegrityError

class FragSampleModelTests(TransactionTestCase):

    def path_and_user(self):
        return 'fragmentation/tests/files/example.mgf', get_user_model().objects.create()

    def test_import_sample(self):
        sample_file_path, u = self.path_and_user()
        with open(sample_file_path, 'r') as fss:
            exp_res = ''.join(fss.readlines())
        with open(sample_file_path, 'rb') as fss:
            fs = FragSample.import_sample(fss, u, 'name', 'file_name', energy=0, task=False)
            fs.wait_import_done()
            self.maxDiff = None
            self.assertEqual(fs.gen_mgf(energy = 0), exp_res)

    def test_ions_limit(self):
        sample_file_path, u = self.path_and_user()
        FragSample.IONS_LIMIT = 2
        with self.assertRaises(IntegrityError):
            with open(sample_file_path, 'rb') as fss:
                fs = FragSample.import_sample(fss, u, 'name', 'file_name', energy=0)
                fs.wait_import_done()
