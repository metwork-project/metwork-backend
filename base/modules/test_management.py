# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import tempfile
from django.test import TransactionTestCase
from django.contrib.auth import get_user_model
from django.test import override_settings


@override_settings(DATA_FILES_PATH=tempfile.gettempdir())
class BaseTestManagement(TransactionTestCase):
    pass


class TestManagement(BaseTestManagement):
    def get_user(self, email="get@user.com"):
        return get_user_model().objects.create(email=email)
