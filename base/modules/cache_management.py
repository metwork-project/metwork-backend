# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
from django.db import models
from django.core.cache import cache

def model_from_cache(cache_key, item_key):
	model_list = cache.get(cache_key.__name__)
	if model_list == None:
		model_list = {}
	if not item_key in model_list:
		model_list[item_key] = cache_key.objects.get(id=item_key)
		cache.set(cache_key.__name__, model_list)
	return model_list[item_key]

#class CacheManagement(object):



