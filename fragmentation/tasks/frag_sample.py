from __future__ import absolute_import, unicode_literals
from celery import shared_task
from fragmentation.models import FragSample
from django.conf import settings


@shared_task
def set_ions_count(ions_count, frag_sample_id):
    fs = FragSample.objects.get(id=frag_sample_id)
    fs.ions_count = sum(ions_count)
    fs.save()
    return fs.id


@shared_task
def import_ion(frag_sample_id, ion, energy):
    fs = FragSample.objects.get(id=frag_sample_id)
    return fs.import_ion(ion, energy)


@shared_task
def gen_cosine_matrix(frag_sample_id):
    fs = FragSample.objects.get(id=frag_sample_id)
    fs.gen_cosine_matrix()
    return fs.id


@shared_task
def get_molecular_network(frag_sample_id):
    fs = FragSample.objects.get(id=frag_sample_id)
    fs.get_molecular_network()
    return fs.id


@shared_task
def finalize_import(frag_sample_id):
    fs = FragSample.objects.get(id=frag_sample_id)
    fs.finalise_import()
    return fs.id
