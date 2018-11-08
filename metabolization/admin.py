# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib import admin
from .models import ReactProcess

# Register your models here.
@admin.register(ReactProcess)
class ReactProcessAdmin(admin.ModelAdmin):
    readonly_fields = ( 'reaction',  'status_code', 'reactants_smiles', 'products_smiles')
    fields = ('reaction', 'status_code', 'reactants_smiles', 'products_smiles')
    list_display = ('id', 'reaction', 'status_code')
    list_filter = ('status_code', )
    search_fields = ['reaction__name']
    def reactants_smiles(self, instance):
        return '.'.join([ r.smiles() for r in instance.reactants.all()])
    reactants_smiles.short_description = 'Reactants'

    def products_smiles(self, instance):
        return '.'.join([ r.smiles() for r in instance.products.all()])
    products_smiles.short_description = 'Products'
