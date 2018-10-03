# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from base.views.model_auth import ModelAuthViewSet
from metabolization.models import Reaction
from rest_framework import serializers
from rest_framework.decorators import list_route, detail_route
from rest_framework.response import Response
from django.contrib.auth import get_user_model
from rest_framework.parsers import JSONParser
from django.http import JsonResponse
from base.modules import JSONSerializerField, ChemDoodle
from base.views import MoleculeSerializer

class ReactionSerializer(serializers.ModelSerializer):
    user = serializers.PrimaryKeyRelatedField(
        queryset = get_user_model().objects.all()
        # ,allow_null = True
        )

    class Meta:
        model = Reaction
        fields = (
            'name',
            'description',
            'user',
            'user_name',
            'reactants_number',
            'has_no_project',
            'status_code',
            'is_reactor',
            'smarts',
            'chemdoodle_json',
            'chemdoodle_json_error',)
            # )

    chemdoodle_json = JSONSerializerField()

class ReactionViewSet(ModelAuthViewSet):
    queryset = Reaction.objects.all().order_by('name')
    serializer_class = ReactionSerializer

    def get_queryset(self):
        queryset = Reaction.objects.all()
        project_id = self.request.query_params.get('project_id', None)
        if project_id is not None:
                from base.models import SampleAnnotationProject
                selected = self.request.query_params.get('selected', None)
                if selected is not None:
                    if selected == 'true':
                        if project_id != '':
                            queryset = SampleAnnotationProject.objects.get(id=project_id).reactions()
                        else:
                            queryset = Reaction.objects.none()
                    else:
                        if project_id != '':
                            queryset = SampleAnnotationProject.objects.get(id=project_id).reactions_not_selected()
        else:
            filter = self.request.query_params.get('filter', None)
            if filter == 'not_obsolete':
                queryset = queryset.filter(status_code__lt = Reaction.status.OBSOLETE)
        return queryset.order_by('name')

    def create(self, request, *args, **kwargs):
        request.data['user'] = request.user.id
        return super().create(request, *args, **kwargs)

    def update(self, request, *args, **kwargs):
        if 'user' in request.data:
            del request.data['user']
        return super().update(request, *args, **kwargs)

    @list_route(methods=['post'])
    def uploadfile(self, request):
        req_data =  request.data
        try:
            fs = Reaction.import_file(
                file_object = req_data['file_data'],
                user = request.user,
                name = req_data['name'],
                description = req_data['description'])
            return Response({'status': 'ok'})
        except:
            return Response({'error': 'unkown error while upploading file'})

    @detail_route(methods=['get'])
    def get_image(self, request, pk=None):
        reaction = self.get_object()
        return Response({'image': reaction.get_image()})

    @detail_route(methods=['patch'])
    def load_smarts(self, request, pk=None):
        try:
            reaction = self.get_object()
            serializer = self.serializer_class(reaction)
            data = JSONParser().parse(request)
            reaction.load_smarts(data['smarts'])
            return JsonResponse({'success': reaction.chemdoodle_json})
        except:
            return JsonResponse({'error': 'error import smarts'})

    @detail_route(methods=['post'])
    def run_reaction(self, request, pk=None):
        data = JSONParser().parse(request)
        chemdoodle_json = data['reactants']['chemdoodle_json']
        if 'm' in chemdoodle_json:
            cd  = ChemDoodle()
            reactants = [ cd.json_to_mol(mol_json)
                for mol_json in chemdoodle_json['m'] ]
            r = self.get_object()
            rp = r.run_reaction(reactants)
            response = {
                'reactants': [ r.chemdoodle_json for r in rp.reactants.all() ],
                'products': [ p.chemdoodle_json for p in rp.products.all() ] }
            return JsonResponse(response)
        else:
            return JsonResponse({'error': 'test'})
