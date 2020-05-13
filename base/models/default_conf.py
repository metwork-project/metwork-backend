# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models


class DefaultConf(models.Model):

    project_class_name = models.CharField(max_length=255, default="")
    app_name = models.CharField(max_length=255, default="")
    conf_class_name = models.CharField(max_length=255, default="")
    conf_default_id = models.IntegerField(default=0)
