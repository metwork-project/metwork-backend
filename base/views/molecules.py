# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from rest_framework import serializers
from base.models import Molecule
from base.views.model_auth import ModelAuthViewSet
from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import csrf_exempt
from django.http import FileResponse

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

@require_http_methods(["GET"])
@csrf_exempt
def download_all_molecules(request):
	file_path = 'all_molecules.csv'
	query_set = Molecule.objects.all()
	Molecule.gen_molecules(file_path, query_set)
	return FileResponse(open(file_path, 'rb'))
