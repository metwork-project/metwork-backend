# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from decimal import *
from base.models import Molecule
from fragmentation.models import FragSample, FragCompareConf, FragSimConf
from fragmentation.modules import FragCompare, FragSim
from django.contrib.auth import get_user_model
from fragmentation.tests.test_frag_common import FragCommonTests

class FragCompareModelTests(TransactionTestCase):

	TEST_FILES_PATH = 'fragmentation/tests/files/'

	@classmethod
	def setUp(cls):
		cls.user = get_user_model().objects.create()

	def new_frag_compare(self):
		return FragCompare(\
			FragCompareConf.objects.create(ppm_tolerance = 5))


	def test_cosine_2_mgf(self):
		#user = get_user_model().objects.create()
		fc = self.new_frag_compare()
		# import fragsamples
		fm = []
		for fn in ['cosine_1.mgf', 'cosine_2.mgf']:
			with open( self.TEST_FILES_PATH + fn , 'rb' ) as fr:
				fs = FragSample.import_sample(fr, FragCompareModelTests.user, 'name', 'file_name', energy=0)
				fs.wait_import_done()
				fm.append(fs.fragmolsample_set.first())
		fmc = fc.compare_frag_mols(fm)
		self.assertNotEqual(fmc, None)

		#self.assertTrue(fmtch.match)
		#print fmc.cosine, fmc.num_frag_match
		#self.assertEqual(fmc.cosine , Decimal('0.863593'), fmc.cosine)
		#self.assertTrue(fmc.cosine >= Decimal('0.3'), fmc.cosine)


# Bug during test if 2 separate tests ....
	def test_cosine_mgf_sim(self):
		#user = get_user_model().objects.create()
		fc = self.new_frag_compare()
		with open( self.TEST_FILES_PATH + 'cosine_4.mgf' , 'rb' ) as fr:
			fsample = FragSample.import_sample(fr, FragCompareModelTests.user, 'name', 'file_name', energy=0)
			fsample.wait_import_done()

		sm = 'N=C(N)NCC/C=C1/N=C(O)N(/C=C/c2ccc(O)cc2)C1=O'
		m = Molecule.load_from_smiles(sm)
		fsim = FragCommonTests.new_frag_sim()
		fms = fsim.frag_molecule(m)
		fmc = fc.compare_frag_mols([fsample.fragmolsample_set.first(), fms])
		#print fmc.cosine, fmc.num_frag_match
		self.assertTrue(fmc.cosine >= Decimal('0.2'), fmc.cosine)








