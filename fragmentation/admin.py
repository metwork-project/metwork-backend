# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib import admin
from .models import FragSample, FragMolSim

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
