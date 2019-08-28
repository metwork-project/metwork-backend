# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import json
from base.views.model_auth import ModelAuthViewSet, IsOwnerOrPublic
from base.modules import TagViewMethods
from fragmentation.models import FragSample
from rest_framework import serializers
from rest_framework.decorators import list_route, detail_route
from rest_framework.response import Response
from django.db import IntegrityError
from django.http import JsonResponse

class FragSampleSerializer(serializers.ModelSerializer):

    class Meta:
        model = FragSample
        fields = (
            'name',
            'file_name',
            'ion_charge',
            'description',
            'tags_list',
            'ions_count',
            'ions_total',
            'annotations_count',
            'has_no_project',
            'status_code')

class FragSampleViewSet(ModelAuthViewSet, TagViewMethods):
    serializer_class = FragSampleSerializer
    queryset = FragSample.objects.all()
    permission_classes = (IsOwnerOrPublic, )

    def get_queryset(self):
        if self.action == 'list':
            return FragSample.objects.filter(user=self.request.user).order_by('-id')
        else:
            return FragSample.objects.all() 

    @list_route(methods=['post'])
    def uploadfile(self, request):
        req_data =  request.data
        try:
            fs = FragSample.import_sample(
                file_object = req_data['file_data'],
                user = request.user,
                file_name = req_data['file_name'],
                name = req_data['name'],
                description = req_data['description'],
                ion_charge = req_data['ion_charge'],
                task=True)
            return Response({'status': 'ok'})
        except IntegrityError as e:
            return Response({'error': str(e)})


    @detail_route(methods=['post'])
    def add_annotation(self, request, pk=None):
        fs = self.get_object()
        req_data =  request.data
        ion_id = req_data['ion_id']
        smiles = req_data['smiles']
        db_source = req_data['db_source']
        db_id = req_data['db_id']
        fa = fs.add_annotation(
            ion_id = ion_id,
            smiles = smiles,
            db_source = db_source,
            db_id = db_id)
        return Response({'status': 'ok'})

    @detail_route(methods=['post'])
    def uploadfile_annotation(self, request, pk=None):
        fs = self.get_object()
        req_data =  request.data
        file_object = req_data['file_data']
        file_format = req_data['file_format']
        try:
            return Response(fs.import_annotation_file(file_object, file_format))
        except:
            return Response({'error': 'unkown error while upploading file'})

    @detail_route(methods=['get'])
    def molecular_network(self, request, pk=None):
        frag_sample = self.get_object()
        data = frag_sample.molecular_network()
        return JsonResponse(data, safe=False)

    @detail_route(methods=['get'])
    def download_mgf(self, request, pk=None):
        fs = self.get_object()

        # Generate .mgf file from database if not exist
        fs.gen_mgf_file()
        
        fileAddress = os.path.join(fs.item_path(), fs.file_name)
        json_data = json.dumps({
            'data': open(fileAddress, 'r').read() })
        return JsonResponse(
            json_data,
            safe=False)