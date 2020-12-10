# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import tempfile
import shutil
from pathlib import Path
from django.test import TransactionTestCase
from django.contrib.auth import get_user_model
from django.test import override_settings

temp_files_path = tempfile.mkdtemp()


@override_settings(DATA_FILES_PATH=temp_files_path)
class BaseTestManagement(TransactionTestCase):
    def _post_teardown(self):
        super()._post_teardown()
        shutil.rmtree(temp_files_path)
        Path(temp_files_path).mkdir()


class TestManagement(BaseTestManagement):
    def get_user(self, email="get@user.com"):
        return get_user_model().objects.create(email=email)
