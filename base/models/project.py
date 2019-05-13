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
					max_length=128,
					default='')
	description = models.CharField(
					max_length=255,
					default='',
					null=True,
					blank=True)
	user = models.ForeignKey(
					settings.AUTH_USER_MODEL,
					on_delete=models.CASCADE,
					db_index=True)
	molecules = models.ManyToManyField(
					Molecule)
	status_code = models.PositiveSmallIntegerField(
					default=0,
					db_index=True)
	public = models.BooleanField(
            default=False)

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

	def user_name(self):
		return self.user.username

	def load_default_conf(self):
		for f in self._meta.local_fields:
			if f.name.endswith('_conf'):
				model = f.related_model
				conf = model()
				query = model.objects.all()
				for conf_field_name in [ lf.name for lf in conf._meta.local_fields if lf.name !='id']:
					query = query.filter(**{ conf_field_name : conf.__getattribute__(conf_field_name) })
				for conf_field in conf._meta.many_to_many:
					query = query.filter( **{ conf_field.name : None } )
				if query.count() > 0:
					conf = query.first()
				else:
					conf.save()
				self.__setattr__( f.name, conf)
		self.save()
		return self

	def update_conf(self,conf_name,params):
		prev_conf = self.__getattribute__(conf_name)
		conf_class = prev_conf.__class__
		query = conf_class.objects.filter(**params)#.filter(**params)
		if query.count() > 0:
			conf = query.first()
		else:
			conf = conf_class.objects.create(**params)
		self.__setattr__(conf_name, conf)
		self.save()
		prev_conf.check_obsolete()
		return self

	def delete(self, *args, **kwargs):
	# When deleting a project, delete conf associated to project
	# if they are not associate with other project.
		prev_confs = [self.__getattribute__(f.name) \
			for f in self._meta.local_fields \
			if hasattr( self.__getattribute__(f.name) , 'IS_CONF') ]
		super(Project, self).delete(*args, **kwargs)
		for conf in prev_confs:
			conf.check_obsolete()

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
			finish_run.apply_async( args= [self.id], queue = settings.CELERY_WEB_QUEUE)
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

	@classmethod
	def achieved_projects(cls):
		return cls.objects\
            .filter(status_code=cls.status.DONE)

	@classmethod
	def running_projects(cls):
		return cls.objects\
            .filter(
				models.Q(status_code=cls.status.QUEUE) \
				| models.Q(status_code=cls.status.RUNNING))
