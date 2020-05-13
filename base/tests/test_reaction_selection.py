# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from django.contrib.auth import get_user_model
from base.models import Molecule, SampleAnnotationProject
from metabolization.models import Reaction
from fragmentation.models import FragSample


class ReactionSelectionTests(TransactionTestCase):
    def test_tag_selection(self):
        # To be Done
        return self.assertTrue(True)
