# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals

from django.db import models
from fragmentation.models import FragMol, FragCompareConf
from decimal import *
from django.core.exceptions import FieldError


class FragMolCompare(models.Model):

    frag_compare_conf = models.ForeignKey(
        FragCompareConf, on_delete=models.PROTECT, default=None, blank=True, null=True
    )
    frag_mols = models.ManyToManyField(FragMol, default=None, blank=True)
    match = models.BooleanField(default=False)
    cosine = models.DecimalField(max_digits=4, decimal_places=3, default=Decimal("0"))
    # energies = "[fragmolId]:[energy],[fragmolId]:[energy]"
    energies = models.CharField(max_length=32, default="")
    num_frag_match = models.IntegerField(default=0)

    def save(self, *args, **kwargs):
        if self.pk and self.frag_mols.count() != 2:
            raise FieldError
        super(FragMolCompare, self).save()
        return self
