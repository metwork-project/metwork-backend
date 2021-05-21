# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from base.models import BaseModel


class Tag(BaseModel):

    name = models.CharField(max_length=64, default="", unique=True)
