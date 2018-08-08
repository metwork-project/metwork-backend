# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from base.models import SampleAnnotationProject
from base.views.model_auth import ModelAuthViewSet
from rest_framework import serializers
from rest_framework.decorators import list_route, detail_route
from rest_framework.response import Response
from django.http import FileResponse, HttpResponse
from rest_framework.parsers import JSONParser
from wsgiref.util import FileWrapper
from django.http import JsonResponse
import os

class ProjectSerializer(serializers.ModelSerializer):
	frag_sample = serializers.PrimaryKeyRelatedField(read_only=True)

	class Meta:
		model = SampleAnnotationProject
		fields = (
			'name',
			'description',
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
	serializer_class = ProjectSerializer
	queryset = SampleAnnotationProject.objects.all()


	def get_queryset(self):
		return SampleAnnotationProject.objects.filter(user=self.request.user).order_by('-id')

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

	@detail_route(methods=['get'])
	def reactions(self, request, pk=None):
		from metabolization.views import ReactionSerializer
		project = self.get_object()
		return Response({'first_test': 'ok'})
		#return Response(ReactionSerializer(project.reactions_conf.reactions.all()).data)

	@detail_route(methods=['patch'])
	def toggle_item(self, request, pk=None):
		project = self.get_object()
		data = JSONParser().parse(request)
		item_id = data['id']
		field = data['field']
		project.toggle_item(field, item_id)
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
	def download_all_molecules(self, request, pk=None):
		project = self.get_object()
		fileAddress = project.item_path() + '/all_molecules.csv'
		if not os.path.isfile(fileAddress):
			project.gen_all_molecules()
		return FileResponse(open(fileAddress, 'rb'))

	@detail_route(methods=['get'])
	def download_annotations(self, request, pk=None):
		project = self.get_object()
		fileAddress = project.item_path() + '/metwork_annotations.csv'
		if not os.path.isfile(fileAddress):
			project.gen_annotations()
		return FileResponse(open(fileAddress, 'rb'))

	@detail_route(methods=['get'])
	def download_annotations_details(self, request, pk=None):
		project = self.get_object()
		fileAddress = project.item_path() + '/metwork_annotations_details.csv'
		if not os.path.isfile(fileAddress):
			project.gen_annotations_details()
		return FileResponse(open(fileAddress, 'rb'))

	@detail_route(methods=['get'])
	def download_metexplore(self, request, pk=None):
		project = self.get_object()
		fileAddress = project.item_path() + '/metexplore.json'
		if not os.path.isfile(fileAddress):
			project.gen_metexplore()
		return FileResponse(open(fileAddress, 'rb'))

	@detail_route(methods=['get'])
	def cytoscapejs_data(self, request, pk=None):
		project = self.get_object()
		data = project.cytoscapejs_data()
		return JsonResponse(data, safe=False)
