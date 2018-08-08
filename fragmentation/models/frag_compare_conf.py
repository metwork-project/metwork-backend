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

    filter_min_intensity = models.FloatField(default=0.0)
    filter_parent_filter_tolerance = models.FloatField(default=0.0)
    filter_matched_peaks_window = models.FloatField(default=0.0)
    filter_min_matched_peaks_search = models.PositiveSmallIntegerField(default=0)

    cosine_mz_tolerance = models.FloatField(default=0.02)
    cosine_min_matched_peaks = models.PositiveSmallIntegerField(default=6)
    cosine_threshold = models.FloatField(default=0.2)
    ppm_tolerance = models.DecimalField(max_digits=6, decimal_places=4, default=5)
