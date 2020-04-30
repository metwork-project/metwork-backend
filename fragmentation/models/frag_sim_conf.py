# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os.path
from django.db import models
from django.conf import settings
from base.modules import ConfManagement

class FragSimConf(ConfManagement, models.Model):

    class JSONAPIMeta:
        resource_name = "fragsimconf"

    PARAM_FILES = ['param', 'conf']
    CFM_ID_FOLDER = settings.CFM_ID_PATH

    threshold = models.DecimalField(max_digits=7, decimal_places=6, default=0.003)
    param_path = models.CharField(max_length=255, default='param/param_output0.log')
    # param_path = param/param_output0.log
    conf_path = models.CharField(max_length=255, default='conf/param_config.txt')
    # conf_path = conf/param_config.txt
    #molecules = models.ManyToManyField(Molecule)

    def file_path(self, file_name):
        return FragSimConf.CFM_ID_FOLDER + getattr(self, file_name + '_path')

    def file_exist(self, file_name):
        return os.path.isfile(self.file_path(file_name))

    def get_fragmol(self, mol):
        return self.fragmol_set.filter(molecule=mol)

    def frag_count(self, met_run):
        return sum(self.has_frag_set(m) for m in met_run.molecules_all())

