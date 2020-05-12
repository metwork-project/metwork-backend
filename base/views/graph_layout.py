# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from rest_framework import serializers
from base.models import GraphLayout
from base.views.model_auth import ModelAuthViewSet, IsOwnerOrPublic


class GraphLayoutSerializer(serializers.ModelSerializer):
    class Meta:
        model = GraphLayout
        # fields = ("name",)


class GraphLayoutViewSet(ModelAuthViewSet):
    # permission_classes = (IsOwnerOrPublic,)
    serializer_class = GraphLayoutSerializer
    queryset = GraphLayout.objects.all()
