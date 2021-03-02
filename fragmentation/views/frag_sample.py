# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import json
from base.views.model_auth import ModelAuthViewSet, IsOwnerOrPublic
from base.modules import TagViewMethods
from fragmentation.models import FragSample
from fragmentation.utils import AnnotationStatus
from rest_framework import serializers
from rest_framework.decorators import action
from rest_framework.response import Response
from django.db import IntegrityError
from django.http import JsonResponse


class FragSampleSerializer(serializers.ModelSerializer):
    class Meta:
        model = FragSample
        fields = (
            "name",
            "file_name",
            "ion_charge",
            "description",
            "tags_list",
            "ions_count",
            "ions_total",
            "annotations_count",
            "has_no_project",
            "status_code",
        )


class FragSampleViewSet(ModelAuthViewSet, TagViewMethods):
    serializer_class = FragSampleSerializer
    queryset = FragSample.objects.all()
    permission_classes = (IsOwnerOrPublic,)

    def get_queryset(self):
        if self.action == "list":
            return FragSample.objects.filter(user=self.request.user).order_by("-id")
        else:
            return FragSample.objects.all()

    @action(detail=False, methods=["post"])
    def uploadfile(self, request):
        req_data = request.data
        try:
            fs = FragSample.import_sample(
                file_object=req_data["file_data"],
                user=request.user,
                file_name=req_data["file_name"],
                name=req_data["name"],
                description=req_data["description"],
                ion_charge=req_data["ion_charge"],
                task=True,
            )
            return Response({"status": "ok"})
        except IntegrityError as e:
            return Response({"error": str(e)})

    @action(detail=True, methods=["post"])
    def add_annotation(self, request, pk=None):
        fs = self.get_object()
        req_data = request.data
        data = {
            key: req_data[key]
            for key in ("ion_id", "smiles", "name", "status_id", "db_source", "db_id")
        }
        fs.add_annotation_from_smiles(**data)
        return Response({"status": "ok"})

    @action(detail=True, methods=["post"])
    def uploadfile_annotation(self, request, pk=None):
        fs = self.get_object()
        req_data = request.data
        file_object = req_data["file_data"]
        file_format = req_data["file_format"]
        try:
            return Response(fs.import_annotation_file(file_object, file_format))
        except:
            return Response({"error": "unkown error while upploading file"})

    @action(detail=True, methods=["get"])
    def molecular_network(self, request, pk=None, force=False):
        frag_sample = self.get_object()
        data = frag_sample.get_molecular_network(force=force, task=True)
        return JsonResponse(data, safe=False)

    @action(detail=True, methods=["get"])
    def download_file(self, request, pk=None):
        fs = self.get_object()
        file_type = request.query_params["file"]
        if file_type == "mgf":
            fs.gen_mgf_file()
            file_name = fs.file_name
        elif file_type == "annotations":
            fs.gen_annotations_file()
            file_name = "annotations.csv"
        file_path = os.path.join(fs.item_path(), file_name)
        json_data = json.dumps({"data": open(file_path, "r").read()})
        return JsonResponse(json_data, safe=False)
