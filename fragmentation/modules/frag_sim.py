# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import subprocess, re
from decimal import *
from fragmentation.models import FragSimConf, FragMol, FragMolSim, FragMolPeak
from django.db.models import Q

class FragSim:

	def __init__(self, frag_sim_conf):
		self.frag_sim_conf = frag_sim_conf

	def get_frag_mol(self, molecule):
		#return FragMol.objects.instance_of(FragMolSim).filter(\
		fm_search = FragMolSim.objects.filter(
					frag_sim_conf = self.frag_sim_conf, 
					molecule = molecule)
		if fm_search.count() > 0:
			return fm_search.first()
		else:
			return False

	def create_frag_mol_sim(self, molecule):
		fms = FragMolSim.objects.create(
				molecule = molecule, 
				frag_sim_conf = self.frag_sim_conf,
				mass = molecule.mass_exact())
		fms.update_hashes()
		return fms

	def frag_molecule(self, molecule):
		get_fm = self.get_frag_mol(molecule)
		if get_fm:
			return get_fm.wait_run_end()
		else:
			fms = self.create_frag_mol_sim(molecule)
			fms.status_code = FragMolSim.status.RUNNING
			fms.save()
			run_out = subprocess.check_output([\
					self.frag_sim_conf.CFM_ID_FOLDER + 'cfm-predict', \
					molecule.smiles(), \
					str(self.frag_sim_conf.threshold), \
					self.frag_sim_conf.file_path('param'), \
					self.frag_sim_conf.file_path('conf') \
				]).decode()
			_en, en = '', None
			for l in run_out.split('\n'):
				_en = re.findall('(?<=energy)(\d)',l, re.U)
				if _en != []:
					en = int(_en[0])
				elif l !='':
					getcontext().prec = FragMol.MASS_DECIMALS
					fp = [Decimal(f) for f in l.split(' ')]
					FragMolPeak.objects.create(\
					frag_mol = fms, \
					energy = en, \
					mz = fp[0], \
					intensity = fp[1] \
					)
			fms.status_code = FragMolSim.status.DONE
			fms.save()
			return fms

