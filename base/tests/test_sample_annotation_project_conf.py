# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import *
from metabolization.models import *
from fragmentation.models import *
from django.contrib.auth import get_user_model
from metabolization.modules import ReactionTestManagement

class SampleAnnotationProjectConfModelTests(ReactionTestManagement):

	def test_default_conf(self):
		project_class_name = "SampleAnnotationProject"
		u = get_user_model().objects.create(email = 'user@test.com')
		p = SampleAnnotationProject.objects.create(
			user = u,)
		for app_name, conf_class_name, p_conf_name in [
			("metabolization", "ReactionsConf", "reactions_conf"),
			("fragmentation", "FragSimConf", "frag_sim_conf"),
			("fragmentation", "FragCompareConf", "frag_compare_conf")]:
			# rc_filter = DefaultConf.objects.filter(
			# 	project_class_name = project_class_name,
			# 	app_name = app_name,
			# 	conf_class_name = conf_class_name)
			# self.assertNotEqual( rc_filter.count(), 0)
			self.assertNotEqual( getattr(p , p_conf_name).id, 0)
		self.assertEqual(
			set([ r.id for r in p.reactions_conf.reactions.all() ]),
			set([ r.id for r in Reaction.objects.all() ]) )
		# test keep custom conf while save
		rc = ReactionsConf.objects.create()
		p.reactions_conf = rc
		p.save()
		self.assertEqual(p.reactions_conf, rc)

	def test_status_ready(self):
		initial_name = "initial name"
		u = get_user_model().objects.create(email = 'user@test.com')
		# self.import_file(reaction_name = "methylation", user = u)
		r = self.create_reacts([('methylation', '[N,O:1]>>[*:1]-[#6]')]) ['methylation']
		p = SampleAnnotationProject.objects.create(
			name = initial_name,
			user = u,)
		self.assertEqual(p.status_code, Project.status.INIT)
		sample_file_path = 'fragmentation/tests/files/example.mgf'
		with open(sample_file_path, 'rb') as fss:
			fs = FragSample.import_sample(fss, u, 'name', 'file_name')
			fs.wait_import_done()
		p.update_frag_sample(fs)
		p.save()
		self.assertEqual(p.status_code, Project.status.INIT)
		fs.add_annotation(1, 'CCC')
		p.update_frag_sample(fs)
		p.save()
		self.assertEqual(p.molecules.count(), 1)
		self.assertEqual(p.status_code, Project.status.INIT)
		rc = ReactionsConf.objects.create()
		self.assertEqual(rc.reactions.count(), 0)
		p.reactions_conf = rc
		p.save()
		self.assertEqual(p.status_code, 0)
		rc.reactions.add(Reaction.objects.first())
		p.save()
		self.assertEqual(p.status_code, Project.status.READY)
		p.run()
		self.assertTrue(p.status_code > Project.status.READY)
		other_name = "other"
		p.name = initial_name + " modified"
		p.save()
		self.assertEqual(p.name, initial_name)
		p.status_code = Project.status.DONE
		p.save()
		self.assertEqual(p.status_code, Project.status.DONE)
		p.status_code = Project.status.READY
		p.save()
		self.assertEqual(p.status_code, Project.status.DONE)

	def clone_project(self):
		#Reaction.reactions_update()
		initial_name = "origin name"
		clone_suffix = " COPY"
		u = get_user_model().objects.create(email = 'user@test.com')
		self.import_file(reaction_name = "methylation", user = u)
		p = SampleAnnotationProject.objects.create(
			name = initial_name,
			user = u,)
		sample_file_path = 'fragmentation/tests/files/example.mgf'
		with open(sample_file_path, 'rb') as fss:
			fs = FragSample.import_sample(fss, u, 'name', 'file_name')
			fs.wait_import_done()
		p.update_frag_sample(fs)
		fs.add_annotation(1, 'CCC')
		p.update_frag_sample(fs)
		rc = ReactionsConf.objects.create()
		p.reactions_conf = rc
		rc.reactions.add(Reaction.objects.first())
		p.save()
		pc = p.clone_project()
		self.assertNotEqual(pc, p)
		self.assertEqual(pc.name, p.name + clone_suffix)
		fields = [
			'user',
			'description',
			'depth_total',
			'depth_last_match',
			'reactions_conf',
			'frag_sim_conf',
			'frag_compare_conf',
			'frag_sample',
		]
		self.assertEqual(pc.user, p.user)
		for f in fields:
			self.assertEqual(getattr(pc, f), getattr(p, f))
		def fais(project):
			return { fai.id for fai in project.frag_annotations_init.all() }
		self.assertEqual(fais(pc), fais(p))
