# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from fragmentation.models import FragMol

class FragMolPeak(models.Model):

    class Meta:
        ordering = ('mz', )

    class JSONAPIMeta:
        resource_name = "fragmolpeak"

    def __str__(self):
        return '{0} {1} energy {2}'.format(self.mz, self.energy, self.energy)

    frag_mol = models.ForeignKey(FragMol, default=None, on_delete=models.CASCADE)#, db_index = True)
    energy = models.SmallIntegerField(default=0)
    mz = models.DecimalField(max_digits=26, decimal_places=10)
    intensity = models.DecimalField(max_digits=26, decimal_places=10)
