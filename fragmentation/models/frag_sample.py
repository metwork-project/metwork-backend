# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from decimal import *
import time
from django.db import models, IntegrityError
#from django.contrib.auth.models import User
from django.conf import settings
#from base.models import User
from base.models import Molecule

class FragSample(models.Model):

	class JSONAPIMeta:
		resource_name = "fragsamples"

	class Meta:
		ordering = ('name', )

	user = models.ForeignKey(
				settings.AUTH_USER_MODEL, 
				on_delete=models.PROTECT, 
				default=None)
	name = models.CharField(
					max_length=128, 
					default='')
	file_name = models.CharField(
					max_length=255, 
					default='')
	description = models.CharField(
					max_length=255, 
					default='',
					null= True,
					blank=True)
	ions_total = models.PositiveSmallIntegerField(
					default=0, 
					db_index = True)
	status_code = models.PositiveIntegerField(
					default=0, 
					db_index = True)

	# Limit the number of ions per sample
	IONS_LIMIT = 100000

	class status:
		INIT = 0
		READY = 1
		RUNNING = 2
		DONE = 3
		ERROR = 99

	def __str__(self):
		return self.name

	def has_no_project(self):
		return self.sampleannotationproject_set.count() == 0

	def ions_count(self):
		return self.fragmolsample_set.count()

	def annotations_count(self):
		from fragmentation.models import FragAnnotationDB
		return FragAnnotationDB.objects.filter(frag_mol_sample__frag_sample = self).count()

	def add_annotation(self, ion_id, smiles, db_source='', db_id=''):
	# ===> Add ion_name !
		from fragmentation.models import FragAnnotationDB
		fms = self.fragmolsample_set.get(ion_id = int(ion_id))
		return FragAnnotationDB.objects.create(
			frag_mol_sample = fms,
			molecule = Molecule.load_from_smiles(smiles),
			db_source = db_source,
			db_id = db_id)

	@classmethod
	def import_sample(cls, file_object, user, name='', file_name='', description='', energy=1, task=False):
		from fragmentation.models import FragMolSample
		from fragmentation.tasks import import_sample_task
		
		data = [l.decode('utf-8') for l in file_object.readlines()]
		#data = file_object.readlines()
		error_log = []
		total_ions = data.count('BEGIN IONS\n')
		
		if total_ions > FragSample.IONS_LIMIT:
			raise IntegrityError(
				'{0} ions max authorized, {1} in the sample.'.format(FragSample.IONS_LIMIT, total_ions))
		fs = FragSample.objects.create(
			user = user, 
			name = name, 
			file_name = file_name, 
			description = description, 
			status_code = 1,
			ions_total = total_ions)
		fs.status_code = 2
		fs.save()
		if task:
			import_sample_task.apply_async(args = [fs.id, data, energy], queue = settings.CELERY_TASK_DEFAULT_QUEUE)
		else:
			fs.import_sample_(data, energy)
		return fs

	def import_sample_(self, data, energy):
		from fragmentation.models import FragMolSample, FragMolAttribute, FragMolPeak
		import re
		error_log = []

		for l in data:
			l = l.replace('\n','')
			#try:
			av = re.split('=', l)
			peak = re.match(u'([\d]*\.+[\d]*) ([\d]*\.+[\d]*)', l, re.U)

			# Create sample mol if begin line
			if l == 'BEGIN IONS':

				fsm = FragMolSample.objects.create(frag_sample = self)
				p = 1

			# Add attribute
			elif len(av) == 2:
				if av[0] == 'PEPMASS':
					fsm.mass = Decimal(av[1])
					fsm.save()
				elif av[0] == 'SCANS': 
					fsm.ion_id = int(av[1])
					fsm.save()
				else:
					FragMolAttribute.objects.create(
						frag_mol = fsm,\
						title = av[0],\
						value = av[1],\
						position = p)
				p +=1
			if peak:
				FragMolPeak.objects.create(
						frag_mol = fsm,
						energy = energy,
						mz = Decimal(peak.group(1)),
						intensity = Decimal(peak.group(2)))
			if l == 'END IONS':
				p = 1
			#except:
			#	error_log.append(l)
		if len(error_log) > 0 : print ('ERROR LOG', error_log )
		self.status_code = 3
		self.save()

	def wait_import_done(self, timeout=360):
		begin = time.time()
		while self.status_code == FragSample.status.RUNNING:
			time.sleep(0.5)
			if (time.time() - begin) > timeout:
				print ('\n#### close due to timeout #####\n')
				return self
			else:
				self.refresh_from_db()
		#time.sleep(20)
		return self

	def import_annotation_file(self, file_object, file_format='default'):
		from fragmentation.models import FragMolSample, FragAnnotationDB
		fls = [l.decode('utf-8') for l in file_object.readlines()]
		#fls = file_object.readlines()
		for fl in fls[1:]:
			if file_format == 'default':
				ion_id, name, smiles, db_source, db_id = fl.split("\n")[0].split(",")
			elif file_format == 'GNPS':
				data = fl.split("\t")
				ion_id, name, smiles, db_source, db_id = data[0], data[4], data[26], 'GNPS : ' + data[5] + ', ' + data[6], data[2]

			if int(ion_id) > 0:
				m = Molecule.load_from_smiles(smiles)
				fms = self.fragmolsample_set.get(
						ion_id = ion_id)

				fa = FragAnnotationDB.objects.create(
					frag_mol_sample = fms,
					molecule = m,
					name = name,
					db_source = db_source,
					db_id = db_id)
		return {'success': 'Annotations successfully imported'}

	def gen_mgf(self, energy = 2, decimal = 6):
		from fragmentation.models import FragMolSample
		res = '\n'.join([\
					fm.gen_mgf(energy) \
					for fm in FragMolSample.objects.filter(frag_sample = self).order_by("ion_id") ])
		return res
