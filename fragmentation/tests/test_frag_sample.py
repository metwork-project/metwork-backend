# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from pathlib import Path
from django.test import TransactionTestCase
from django.contrib.auth.models import User
from fragmentation.models import FragSample
from django.contrib.auth import get_user_model
from django.db import IntegrityError


class FragSampleTestUtils:

    user = None
    sample_file_path = "fragmentation/tests/files/example.mgf"

    def get_user(self):
        if not self.user:
            self.user = get_user_model().objects.create()

    def create_frag_sample(self, sample_file_path=None):
        if not sample_file_path:
            sample_file_path = self.sample_file_path
        self.get_user()
        with open(sample_file_path, "rb") as fss:
            fs = FragSample.import_sample(
                fss, self.user, "name", "file_name", task=False
            )
            fs.wait_import_done()
            return fs


class FragSampleModelTests(TransactionTestCase, FragSampleTestUtils):
    def test_import_sample(self):
        exp_res = Path(self.sample_file_path).read_text()
        fs = self.create_frag_sample()
        self.assertEqual(fs.gen_mgf(energy=1), exp_res)

    def test_ions_limit(self):
        FragSample.conf["IONS_LIMIT"] = "2"
        with self.assertRaises(IntegrityError):
            self.create_frag_sample()
        FragSample.IONS_LIMIT = 10000

    def test_get_molecular_network(self):
        fs = self.create_frag_sample()
        assert fs.get_molecular_network()
