# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from polymorphic.models import PolymorphicModel
from base.models import Molecule
from fragmentation.models import FragMolSample, FragMolCompare

class FragAnnotation(PolymorphicModel):

    frag_mol_sample = models.ForeignKey(
            FragMolSample,
            on_delete=models.CASCADE,
            default=None)
    molecule = models.ForeignKey(
            Molecule,
            on_delete=models.PROTECT)

    def ion_id(self):
        return self.frag_mol_sample.ion_id

    def smiles(self):
        return self.molecule.smiles()

    def chemdoodle_json(self):
        return self.molecule.chemdoodle_json

class FragAnnotationDB(FragAnnotation):

    class JSONAPIMeta:
        resource_name = "frag-annotations"

    name = models.CharField(
            max_length=256,
            default='')
    db_source = models.CharField(
            max_length=256,
            default='unkown')
    db_id = models.CharField(
            max_length=128,
            default='')

    def has_no_project(self):
        return self.sampleannotationproject_set.count() == 0

    def delete(self, *args, **kwargs):
        if self.has_no_project():
            super(FragAnnotationDB, self).delete(*args, **kwargs)

class FragAnnotationCompare(FragAnnotation):

    project = models.ForeignKey(
            'base.Project',
            on_delete=models.CASCADE,
            db_index = True,
            default=None)
    frag_mol_compare = models.ForeignKey(
            FragMolCompare,
            on_delete=models.PROTECT,
            default=None,
            null=True)
