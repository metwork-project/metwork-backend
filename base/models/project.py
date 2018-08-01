# -*- coding: utf-8 -*-

# Base class for projects, use polymorhpic models so more specfic projects
# as SampleAnnotationsProject are based on this class

from __future__ import unicode_literals
import time
from django.db import models
from django.apps import apps
from polymorphic.models import PolymorphicModel
from base.modules import FileManagement
from django.conf import settings
from base.models import Molecule
from django.core.cache import cache

class Project(FileManagement, PolymorphicModel):

	name = models.CharField(
					max_length=64, 
					default='')
	description = models.CharField(
					max_length=255, 
					default='',
					null= True,
					blank=True)
	user = models.ForeignKey(
					settings.AUTH_USER_MODEL, 
					on_delete=models.CASCADE, 
					db_index = True)
	molecules = models.ManyToManyField(
					Molecule)
	status_code = models.PositiveSmallIntegerField(
					default=0, 
					db_index = True)

	class JSONAPIMeta:
		resource_name = "projects"

	class status:
		INIT = 0
		READY = 1
		QUEUE = 2
		RUNNING = 3
		DONE = 4
		ERROR = 99

	def __str__(self):
		return self.name

	def load_default_conf(self):
	# give the default confs for *_conf filed
	# create default if not exist
	# ===> To be simplified !
		from base.models import DefaultConf
		project_class_name = self.__class__.__name__
		for f in self._meta.local_fields:
			f_name = f.name
			if f_name.endswith('_conf'):
				#default_conf = None # ??
				conf_default_id = 0
				model = f.related_model
				app_name = model._meta.app_label
				conf_class_name = model.__name__
				conf_class = apps.get_model(app_name, conf_class_name)
				id_lookup =  DefaultConf.objects\
							.filter( 
								project_class_name = project_class_name,
								app_name = app_name,
								conf_class_name = conf_class_name)
				id_lookup_count = id_lookup.count()
				if id_lookup_count > 0:
					conf_default_id = id_lookup.first().conf_default_id
				if conf_default_id == 0:
					if conf_class.objects.count() == 0:
						default_conf = conf_class.objects.create()
					else:
						default_conf = conf_class.objects.first()
					conf_default_id = default_conf.id
				else:
					default_conf = conf_class.objects.get(id = conf_default_id)
				setattr(self, f_name, default_conf)
				if id_lookup_count == 0: 
					default_conf = DefaultConf(
							project_class_name = project_class_name,
							app_name = app_name,
							conf_class_name = conf_class_name)
				default_conf.conf_default_id = conf_default_id
				default_conf.save()
		self.save()
		return self

	def delete(self, *args, **kwargs):
	# When deleting a project, delete conf associated to project 
	# if they are not associate with other project.
	# ===> To add : not delete if default conf
		confs = [(f.related_model, self.__getattribute__(f.attname)) \
			for f in self._meta.local_fields if f.name.endswith('_conf')]
		super(Project, self).delete(*args, **kwargs)
		for cl, conf_id in confs:
			conf = cl.objects.get(id=conf_id)
			for ro in conf._meta.related_objects:
				accessor = ro.get_accessor_name()
				if accessor.endswith('project_set') \
					and conf.__getattribute__(accessor).count() == 0:
					try:
						conf.delete()
					except:
						pass

	def process_count_key(self):
		return 'project_process_' + str(self.id)

	def add_process(self, count=1):
		cache.incr(self.process_count_key(), count)
		return self

	def close_process(self):
		from base.tasks.project import finish_run
		key = self.process_count_key()
		cache.decr(key)
		if cache.get(key) == 0:
			finish_run.apply_async( args= [self.id], queue = settings.CELERY_TASK_DEFAULT_QUEUE)
			cache.delete(key)
		return self

	def finish_run(self):
		self.status_code = Project.status.DONE
		self.save()

	def wait_run_end(self, timeout = 360):
		# Used for test only
		# timeout is in seconds
		begin = time.time()
		while self.status_code in [Project.status.QUEUE, Project.status.RUNNING]:
			time.sleep(0.5)
			if (time.time() - begin) > timeout:
			# ===> Manage Errors here
				print ('\n#### close due to timeout #####\n')
				return self
			else:
				self.refresh_from_db()
		return self

