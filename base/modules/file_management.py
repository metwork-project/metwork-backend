# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.db import models
import os, shutil
from django.conf import settings


class FileManagement(object):
    def model_path(self):
        return "/".join([settings.DATA_FILES_PATH, self.JSONAPIMeta.resource_name])

    def item_path(self):
        return "/".join([self.model_path(), str(self.pk)])

    def gen_item(self):
        for d in [self.model_path(), self.item_path()]:
            if not os.path.isdir(d):
                os.mkdir(d)
        return self

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        if settings.EDIT_FILES:
            self.gen_item()
        return self

    def delete(self, *args, **kwargs):
        d = self.item_path()
        if os.path.isdir(d):
            # try: ???
            shutil.rmtree(d)
        super().delete(*args, **kwargs)
