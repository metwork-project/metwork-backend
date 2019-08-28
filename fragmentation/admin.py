# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib import admin
from .models import FragSample, FragMolSim, FragMolCompare

@admin.register(FragSample)
class FragSampleAdmin(admin.ModelAdmin):
    readonly_fields = ( 'ions_total', 'projects_count')
    fields = ('name', 'user', 'status_code', 'ions_total', 'projects_count')
    list_display = ('name', 'user', 'status_code', 'ions_total', 'projects_count')
    list_filter = ('status_code', 'user')

    def projects_count(self, instance):
        return instance.sampleannotationproject_set.count()
    projects_count.short_description = 'Projects'

@admin.register(FragMolSim)
class FragMolSimAdmin(admin.ModelAdmin):
    readonly_fields = ( 'smiles',  'status_code')
    fields = ('smiles', 'status_code')
    list_display = ('id', 'smiles', 'status_code')
    list_filter = ('status_code', )

    def smiles(self, instance):
        return instance.molecule.smiles()
    smiles.short_description = 'Smiles'

@admin.register(FragMolCompare)
class FragSampleAdmin(admin.ModelAdmin):
    readonly_fields = ('match', 'cosine', 'str_frag_mols')
    fields = ('match', 'cosine', 'str_frag_mols')
    list_display = ('id', 'str_frag_mols')
    list_filter = ('cosine',)

    def frag_mol_label(self, frag_mol):
        if frag_mol.__class__.__name__ == "FragMolSim":
            return frag_mol.molecule.smiles()
        if frag_mol.__class__.__name__ == "FragMolSample":
            return 'sample: {0}, ion_id: {1}'.format(
                frag_mol.frag_sample.id,
                frag_mol.ion_id)
        return 'None'

    def str_frag_mols(self, instance):
        return ', '.join([self.frag_mol_label(fm) for fm in instance.frag_mols.all()])
    str_frag_mols.short_description = 'Mols'
