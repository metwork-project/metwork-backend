# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from django.contrib.postgres.fields import ArrayField

class Array1DModel(models.Model):

    value = ArrayField(
        models.FloatField(),
        null = True)

class Array2DModel(models.Model):
    value = ArrayField(
        ArrayField(models.FloatField()),
        null = True)
