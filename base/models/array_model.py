# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from django.contrib.postgres.fields import ArrayField
from base.models import BaseModel


class Array1DModel(BaseModel):

    value = ArrayField(models.FloatField(), null=True)


class Array2DModel(BaseModel):
    value = ArrayField(ArrayField(models.FloatField()), null=True)
