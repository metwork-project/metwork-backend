from __future__ import unicode_literals

from django.db import models
from fragmentation.models import FragMol
from django.contrib.postgres.fields import ArrayField
from base.models import BaseModel


class FragMolSpectrum(BaseModel):
    class JSONAPIMeta:
        resource_name = "fragmolspectrum"

    frag_mol = models.ForeignKey(
        FragMol, default=None, on_delete=models.CASCADE
    )  # , db_index = True)
    energy = models.SmallIntegerField(default=0)
    spectrum = ArrayField(ArrayField(models.FloatField()))
