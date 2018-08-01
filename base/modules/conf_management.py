# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.db import models
import os, shutil
from django.conf import settings


class ConfManagement(object):
        
    def delete(self, *args, **kwargs):
        from base.models import DefaultConf
        DefaultConf.objects.filter(
            conf_class_name = self.__class__.__name__,
            conf_default_id = self.id)\
            .delete()
        super(ConfManagement, self).delete(*args, **kwargs)
