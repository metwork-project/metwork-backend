# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from base.views.model_auth import ModelAuthViewSet
from fragmentation.models import FragAnnotationDB
from rest_framework import serializers
from rest_framework.response import Response
from base.modules import JSONSerializerField

class FragAnnotationSerializer(serializers.ModelSerializer):

    class Meta:
        model = FragAnnotationDB
        fields = (
            'ion_id',
            'name',
            'smiles',
            'db_source',
            'db_id',
            'has_no_project',
            'chemdoodle_json',
			'adduct')

    chemdoodle_json = JSONSerializerField()

class FragAnnotationViewSet(ModelAuthViewSet):

	serializer_class = FragAnnotationSerializer
	queryset = FragAnnotationDB.objects.all()

	def get_queryset(self):
		queryset = FragAnnotationDB.objects.filter(
			frag_mol_sample__frag_sample__user=self.request.user)
		if 'frag_sample_id' in self.request.query_params:
			queryset = queryset.filter(frag_mol_sample__frag_sample=self.request.query_params['frag_sample_id'])
		project_id = self.request.query_params.get('project_id', None)
		if project_id is not None:
				from base.models import SampleAnnotationProject
				selected = self.request.query_params.get('selected', None)
				if selected is not None:
					if selected == 'true':
						if project_id != '':
							queryset = SampleAnnotationProject.objects.get(id=project_id).frag_annotations_init.all()
						else:
							queryset = FragAnnotationDB.objects.none()
					else:
						if project_id != '':
							queryset = SampleAnnotationProject.objects.get(id=project_id).frag_annotations_init_not_selected()
		return queryset.order_by('frag_mol_sample__ion_id')
