# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from base.views.model_auth import ModelAuthViewSet
from fragmentation.models import FragCompareConf
from rest_framework import serializers

class FragCompareConfSerializer(serializers.ModelSerializer):

    class Meta:
        model = FragCompareConf
        fields = (
            'filter_min_intensity',
            'filter_parent_filter_tolerance',
            'filter_matched_peaks_window',
            'filter_min_matched_peaks_search',
            'cosine_mz_tolerance',
            'cosine_min_matched_peaks',
            'cosine_threshold' )

class FragCompareConfViewSet(ModelAuthViewSet):
    serializer_class = FragCompareConfSerializer
    queryset = FragCompareConf.objects.all()

    # def get_object(self):
    #     queryset = self.get_object()
    #     project_id = self.request.query_params.get('project_id', None)
    #     print(project_id)
    #     if project_id is not None:
    #         from base.models import SampleAnnotationProject
    #         print(SampleAnnotationProject.objects\
    #                     .get(id = project_id).frag_compare_conf)
    #         queryset = SampleAnnotationProject.objects\
    #                     .filter(id = project_id).frag_compare_conf
    #     return queryset
