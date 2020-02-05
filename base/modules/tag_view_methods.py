# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from rest_framework import viewsets
from rest_framework.parsers import JSONParser
from django.http import JsonResponse
from base.models import  Tag
from rest_framework.decorators import action

class TagViewMethods(viewsets.ModelViewSet):

    @action(detail=True, methods=['patch'])
    def add_tag(self, request, pk=None):
        # try:
            target = self.get_object()
            tag_name = JSONParser().parse(request)
            tag_find = Tag.objects.filter(name=tag_name)
            if tag_find.count() == 1:
                tag = tag_find.first()
            else:
                tag = Tag.objects.create(name=tag_name)
            target.tags.add(tag)
            return JsonResponse({'success': tag_name})
        # except:
        #     return JsonResponse({'error': 'error import smarts'})

    @action(detail=True, methods=['patch'])
    def remove_tag(self, request, pk=None):
        # try:
            target = self.get_object()
            tag_name = JSONParser().parse(request)
            tag_find = target.tags.filter(name=tag_name)
            if tag_find.count() == 1:
                tag = tag_find.first()
                target.tags.remove(tag)
                print(dir(tag))
                if (tag.reaction_tags.count() + tag.fragsample_tags.count()) == 0:
                    tag.delete()
            return JsonResponse({'success': tag_name})
        # except:
        #     return JsonResponse({'error': 'error import smarts'})