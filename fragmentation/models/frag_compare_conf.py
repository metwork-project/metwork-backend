# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models

import os.path
import re
from base.models import Molecule
from base.modules import ConfManagement

class FragCompareConf(ConfManagement, models.Model):

    class JSONAPIMeta:
        resource_name = "fragcompareconf"

    cosine_threshold = models.DecimalField(max_digits=4, decimal_places=3, default=0.18)
    ppm_tolerance = models.DecimalField(max_digits=6, decimal_places=4, default=5)

