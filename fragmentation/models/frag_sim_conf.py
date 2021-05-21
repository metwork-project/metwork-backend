# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os.path
from pathlib import Path
from django.db import models
from django.conf import settings
from base.modules import ConfManagement
from base.models import BaseModel


class FragSimConf(ConfManagement, BaseModel):
    class JSONAPIMeta:
        resource_name = "fragsimconf"

    PARAM_FILES = ["param", "conf"]
    CFM_ID_FOLDER = settings.CFM_ID_PATH
    PARAM_PATH_BY_CHARGE = {
        "positive": "param/param_output0.log",
        "negative": "param/param_output0_neg.log",
    }
    CONF_PATH_BY_CHARGE = {
        "positive": "conf/param_config.txt",
        "negative": "conf/param_config_neg.txt",
    }

    threshold = models.DecimalField(max_digits=7, decimal_places=6, default=0.003)
    param_path = models.CharField(
        max_length=255, default=PARAM_PATH_BY_CHARGE["positive"]
    )
    # param_path = param/param_output0.log
    conf_path = models.CharField(
        max_length=255, default=CONF_PATH_BY_CHARGE["positive"]
    )
    # conf_path = conf/param_config.txt
    # molecules = models.ManyToManyField(Molecule)

    def file_path(self, file_name):
        path = Path(getattr(self, file_name + "_path"))
        if not path.is_absolute():
            path = FragSimConf.CFM_ID_FOLDER / path
        return str(path)

    def file_exist(self, file_name):
        return os.path.isfile(self.file_path(file_name))

    def get_fragmol(self, mol):
        return self.fragmol_set.filter(molecule=mol)

    def frag_count(self, met_run):
        return sum(self.has_frag_set(m) for m in met_run.molecules_all())

    @classmethod
    def params_for_ion_charge(cls, ion_charge):
        return {
            "param_path": cls.PARAM_PATH_BY_CHARGE[ion_charge],
            "conf_path": cls.CONF_PATH_BY_CHARGE[ion_charge],
        }
