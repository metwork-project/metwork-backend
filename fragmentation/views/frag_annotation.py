# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.http import JsonResponse
from rest_framework import serializers
from rest_framework.response import Response
from rest_framework.decorators import action
from base.views.model_auth import IsOwnerOrPublic
from base.views.utils import MetaIdsViewSet
from fragmentation.models import FragAnnotationDB
from base.modules import JSONSerializerField
from base.modules.queryset import FilteredQueryset


class FragAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = FragAnnotationDB
        fields = (
            "ion_id",
            "status_id",
            "mz",
            "name",
            "smiles",
            "db_source",
            "db_id",
            "has_no_project",
            "chemdoodle_json",
            "adduct",
        )

    chemdoodle_json = JSONSerializerField()


class FragAnnotationQueryset(FilteredQueryset):
    def filter_init(self):
        frag_sample_id = self.request.query_params.get("frag_sample_id", None)
        if not (self.project_id or frag_sample_id):
            raise Exception
        if frag_sample_id:
            self.queryset = self.queryset.filter(
                frag_mol_sample__frag_sample=frag_sample_id
            )

    def get_default(self, queryset_):
        return queryset_.get_all_fragsample_annotations()

    def get_all(self, queryset_):
        print("queryset_", queryset_)
        return queryset_.all_annotations_init()

    def get_notselected(self, queryset_):
        return queryset_.frag_annotations_init_not_selected()

    def filter_other(self):
        queryset = self.queryset

        for key, value in self.request.query_params.lists():

            if key == "filter[text]":
                text = value[0]
                if text:
                    queryset = queryset.filter(
                        Q(name__icontains=text) | Q(tags__name__icontains=text)
                    )
            elif key == "filter[my]":
                my = value[0].lower() == "true"
                if my:
                    queryset = queryset.filter(user=self.request.user)
            elif key == "filter[user]" and not my:
                user_text = value[0]
                if user_text:
                    queryset = queryset.filter(user__username__icontains=user_text)

        self.queryset = queryset

    def _filter_status(self, filter_status):
        self.queryset = self.queryset.filter(status_id__in=filter_status)


class FragAnnotationViewSet(MetaIdsViewSet):

    serializer_class = FragAnnotationSerializer
    permission_classes = (IsOwnerOrPublic,)

    def get_queryset(self):
        try:
            queryset = FragAnnotationQueryset(self.request, FragAnnotationDB).queryset
            return queryset.order_by("frag_mol_sample__ion_id")
        except:
            return FragAnnotationDB.objects.filter(id=0)

    @action(detail=True, methods=["get"])
    def get_mgf(self, request, pk=None):
        frag_annot = self.get_object()
        mgf = frag_annot.frag_mol_sample.gen_mgf()
        return JsonResponse({"mgf": mgf}, safe=False)
