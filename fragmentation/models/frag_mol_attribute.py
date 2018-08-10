# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from fragmentation.models import FragMol

class FragMolAttribute(models.Model):

    frag_mol = models.ForeignKey(FragMol, default=None, on_delete=models.CASCADE)
    title = models.CharField(max_length=32, default='')
    value = models.CharField(max_length=128, default='')
    position = models.PositiveSmallIntegerField(default=0)
