# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from rest_framework import serializers
from base.models import Tag
from base.views.model_auth import ModelAuthViewSet

class TagSerializer(serializers.ModelSerializer):
    class Meta:
        model = Tag
        fields = ('name',)

class TagViewSet(ModelAuthViewSet):
    serializer_class = TagSerializer
    queryset = Tag.objects.all()