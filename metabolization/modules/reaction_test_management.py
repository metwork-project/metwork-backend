# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.contrib.auth import get_user_model
from base.modules import TestManagement
from metabolization.models import Reaction, ReactProcess


class ReactionTestManagement(TestManagement):
    def create_reacts(self, reacts, email="create@react.com", activate=True):
        user = self.get_user(email)
        result = {}
        kwargs = {"user": user}
        if activate:
            kwargs["status_code"] = Reaction.status.ACTIVE
        for name, smarts in reacts:
            kwargs.update({"name": name, "smarts": smarts})
            react = Reaction.objects.create(**kwargs)
            result[name] = react
        return result
