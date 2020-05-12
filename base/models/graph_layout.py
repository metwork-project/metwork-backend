# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from django.contrib.postgres.fields import JSONField
from django.conf import settings


class GraphLayout(models.Model):
    class types:
        FRAG_SAMPLE = 1
        METABOLIZATION = 2

    user = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, db_index=True
    )
    graph_type = models.IntegerField()
    data_id = models.IntegerField()
    layout = JSONField()
