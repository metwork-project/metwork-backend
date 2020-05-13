# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.contrib.auth import get_user_model
from base.modules import TestManagement
from metabolization.models import Reaction, ReactProcess


class ReactionTestManagement(TestManagement):
    def create_reacts(self, reacts, email="create@react.com"):
        user = self.get_user(email)
        rd = {}
        for name, smarts in reacts:
            r = Reaction.objects.create(user=user, name=name, smarts=smarts)
            rd[name] = r
        return rd
