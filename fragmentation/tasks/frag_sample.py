from __future__ import absolute_import, unicode_literals
from celery import shared_task
from fragmentation.models import FragSample, FragMolSample, FragMolAttribute, FragMolPeak
import re
from decimal import *
from django.conf import settings


@shared_task
def import_sample_task(frag_sample_id, data, energy):
	fs = FragSample.objects.get(id=frag_sample_id)
	error_log = []

	fs.import_sample_( data, energy, True)

@shared_task
def gen_cosine_matrix_task(frag_sample_id):
	fs = FragSample.objects.get(id=frag_sample_id)
	fs.gen_cosine_matrix()

