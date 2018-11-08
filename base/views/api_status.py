# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from rest_framework import serializers
from rest_framework.permissions import AllowAny
from base.views.model_auth import ModelAuthViewSet
from base.models import APIStatus
from rest_framework import viewsets

class APIStatusSerializer(serializers.ModelSerializer):
    class Meta:
        model = APIStatus
        fields = (
            'available',
            'molecules_count',
            'achieved_projects_count',
            'running_projects_count')

class APIStatusViewSet(viewsets.ModelViewSet):
    serializer_class = APIStatusSerializer
    queryset = APIStatus.objects.all()
    permission_classes = (AllowAny,)

    def get_queryset(self):
        if APIStatus.objects.count() == 0:
            status = APIStatus.objects.create()
        return APIStatus.objects.all()
