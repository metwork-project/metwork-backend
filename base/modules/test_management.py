# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.test import TransactionTestCase
from django.contrib.auth import get_user_model

class TestManagement(TransactionTestCase):

    def get_user(self, email='get@user.com'):
        return get_user_model().objects.create(email=email)
