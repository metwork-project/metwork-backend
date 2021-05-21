# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from metabolization.models import Reaction
from base.modules import ConfManagement
from base.models import BaseModel


class ReactionsConf(ConfManagement, BaseModel):
    class JSONAPIMeta:
        resource_name = "reactconfs"

    METHODS_CHOICES = (
        ("reaction", "Reaction"),
        ("reactor", "Reactor"),
        ("rdkit", "RDKit"),
    )

    reactions = models.ManyToManyField(Reaction)
    method_priority = models.CharField(
        max_length=32, choices=METHODS_CHOICES, default="reaction"
    )  # cls. methods_allowed
