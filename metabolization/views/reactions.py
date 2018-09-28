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
from base.modules import JSONSerializerField

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
        filter = self.request.query_params.get('filter', None)
        if filter == 'is_active':
            queryset = queryset.filter(status_code = Reaction.status.ACTIVE)
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
            return Response({'success': reaction.get_chemdoodle_json()})
        except:
            return Response({'error': 'error import smarts'})

    @detail_route(methods=['get'])
    def run_reaction(self, request, pk=None):
        from base.models import Molecule
        from metabolization.models import ReactProcess
        r = self.get_object()
        reactant_smiles = request.query_params.get('smiles','')
        reactants = [ Molecule.load_from_smiles(sm) for sm in reactant_smiles.split('\n') ]
        rp = ReactProcess.objects.create()
        rp.reaction = r
        rp.reactants.set(reactants)
        rp.save()
        rp.run_reaction()
        return Response({'products': [p.smiles() for p in rp.products.all()]})
