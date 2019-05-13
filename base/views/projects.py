# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from base.models import SampleAnnotationProject
from base.views.model_auth import ModelAuthViewSet, IsOwnerOrPublic
from rest_framework import serializers
from rest_framework.decorators import list_route, detail_route
from rest_framework.response import Response
from django.http import FileResponse, HttpResponse
from rest_framework.parsers import JSONParser
from wsgiref.util import FileWrapper
from django.http import JsonResponse
import os, json

class ProjectSerializer(serializers.ModelSerializer):
    frag_sample = serializers.PrimaryKeyRelatedField(read_only=True)

    class Meta:
        model = SampleAnnotationProject
        fields = (
            'name',
            'description',
            'public',
            'frag_sample',
            'status_code',
            'reaction_ids',
            'REACTIONS_LIMIT',
            'annotation_init_ids',
            'depth_total',
            'depth_last_match',
            'molecules_matching_count',
            'molecules_all_count',
            'frag_compare_conf_id', )

class ProjectViewSet(ModelAuthViewSet):
    permission_classes = (IsOwnerOrPublic, )
    serializer_class = ProjectSerializer


    def get_queryset(self):
        if self.action == 'list':
            return SampleAnnotationProject.objects \
                .filter(user=self.request.user).order_by('-id')
        else:
            return SampleAnnotationProject.objects.all()

    def perform_create(self, serializer):
        #print self.request.user
        serializer.save(user=self.request.user)

    @detail_route(methods=['post'])
    def clone_project(self, request, pk=None):
        project = self.get_object()
        try:
            clone = project.clone_project()
            return Response({'clone_id': clone.id})
        except:
            return Response({'error': 'error while cloning project'})

    @detail_route(methods=['patch'])
    def update_frag_sample(self, request, pk=None):
        from fragmentation.models import FragSample
        project = self.get_object()
        data = JSONParser().parse(request)
        fs = FragSample.objects.get(id = data['frag_sample_id'])
        if fs.user == self.request.user:
            project.update_frag_sample(fs)
        return Response(ProjectSerializer(project).data)

# To delete ???
    @detail_route(methods=['get'])
    def reactions(self, request, pk=None):
        from metabolization.views import ReactionSerializer
        project = self.get_object()
        return Response({'first_test': 'ok'})
        #return Response(ReactionSerializer(project.reactions_conf.reactions.all()).data)

    def change_item(self, self_, request, func):
        project = self_.get_object()
        data = JSONParser().parse(request)
        dataLabel = data['dataLabel']
        if 'item_ids' in data:
            getattr(project, func)(dataLabel, data['item_ids'])
        else:
            getattr(project, func)(dataLabel)
        return Response({'project_id': project.id})

    @detail_route(methods=['patch'])
    def add_items(self, request, pk=None):
        return self.change_item(self, request, 'add_items')

    @detail_route(methods=['patch'])
    def add_all(self, request, pk=None):
        return self.change_item(self, request, 'add_all')

    @detail_route(methods=['patch'])
    def remove_all(self, request, pk=None):
        return self.change_item(self, request, 'remove_all')

    @detail_route(methods=['patch'])
    def remove_item(self, request, pk=None):
        return self.change_item(self, request, 'remove_item')

    @detail_route(methods=['patch'])
    def select_reactions_by_mass(self, request, pk=None):
        project = self.get_object()
        project.select_reactions_by_mass()
        return Response({'project_id': project.id})

    @detail_route(methods=['patch'])
    def update_frag_compare_conf(self, request, pk=None):
        project = self.get_object()
        params = JSONParser().parse(request)
        project.update_conf('frag_compare_conf',params)
        return Response({'frag_compare_conf':project.frag_compare_conf.id})

    @detail_route(methods=['post'])
    def start_run(self, request, pk=None):
        project = self.get_object()
        project.save()
        project.run()
        return Response({'status_code':project.status_code})

    @detail_route(methods=['get'])
    def download_file(self, request, pk=None):
        file = self.request.query_params.get('file', None)
        file_name = self.request.query_params.get('file_name', None)
        methods = {
            'annotations' : 'gen_annotations',
            'annotations_details' : 'gen_annotations_details',
            'metexplore' : 'gen_metexplore',
            'all_molecules' : 'gen_all_molecules',
        }
        project = self.get_object()
        fileAddress = project.item_path() + '/' + file_name
        if not os.path.isfile(fileAddress):
            getattr(project, methods[file])()
        json_data = json.dumps({
            'data': open(fileAddress, 'r').read() })
        return JsonResponse(
            json_data,
            safe=False)

    @detail_route(methods=['get'])
    def metabolization_network(self, request, pk=None):
        project = self.get_object()
        data = project.metabolization_network()
        return JsonResponse(data, safe=False)
