# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.db import models
import os, shutil
from django.conf import settings


class ConfManagement(object):

    def check_obsolete(self):
        if self.project_set.count() == 0:
            self.delete()
