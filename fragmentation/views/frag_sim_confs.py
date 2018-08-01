# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from base.views.model_auth import ModelAuthViewSet
from fragmentation.models import FragSimConf
from fragmentation.serializers import FragSimConfSerializer
from rest_framework.decorators import detail_route
from rest_framework.response import Response
from django.http import HttpResponse
from base import tasks

class FragSimConfViewSet(ModelAuthViewSet):
    serializer_class = FragSimConfSerializer
    queryset = FragSimConf.objects.all()
