# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models


class Tag(models.Model):

    name = models.CharField(max_length=64, default="", unique=True)
