# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from rest_framework import serializers
from base.models import Molecule
from base.views.model_auth import ModelAuthViewSet
from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import csrf_exempt
from rest_framework.decorators import action
from django.http import HttpResponse, JsonResponse
from rest_framework.parsers import JSONParser

class MoleculeSerializer(serializers.ModelSerializer):
    class Meta:
        model = Molecule
        fields = (
            'smiles',
            'chemdoodle_json')

class MoleculeViewSet(ModelAuthViewSet):
    serializer_class = MoleculeSerializer
    queryset = Molecule.objects.all()
    def get_queryset(self):
        met_run_in = self.request.query_params.get('met_run_in', None)
        if met_run_in is not None:
            return Molecule.met_run_in(met_run_in)
        met_run_out = self.request.query_params.get('met_run_out', None)
        if met_run_out is not None:
            return Molecule.met_run_out(met_run_out)
        return self.queryset

    @action(detail=False, methods=['patch'])
    def load_smiles(self, request, pk=None):
        data = JSONParser().parse(request)
        mols = Molecule.load_from_smiles(data['smiles'].replace('\n','.'))
        if mols != False:
            return JsonResponse({'success': [mols.chemdoodle_json]})
        else:
            return JsonResponse({'error': 'unable to load smiles'})
