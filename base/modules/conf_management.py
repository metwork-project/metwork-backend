# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.db import models
import os, shutil
from django.conf import settings
from django.db.models.fields import Field


class ConfManagement(object):
    def check_obsolete(self):
        if self.sampleannotationproject_set.count() == 0:
            self.delete()

    def get_as_dict(self):
        return {
            f.name: getattr(self, f.name)
            for f in getattr(self, "_meta").get_fields()
            if isinstance(f, Field)
        }
