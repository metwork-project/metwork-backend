# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import Molecule
from decimal import *
from fragmentation.tests.test_frag_common import FragCommonTests

import hashlib

class FragMolModelTests(TransactionTestCase):

  ## FragSet

    def new_frag_mol_sim(self, smiles = 'N=C([NH3+])', adduct='M+'):
        fs = FragCommonTests.new_frag_sim()
        m = Molecule.load_from_smiles(smiles)
        return fs.create_frag_mol_sim(molecule=m, adduct=adduct), fs

    def test_create_frag_set(self):
        fms, fs = self.new_frag_mol_sim()
        for fn in fs.frag_sim_conf.PARAM_FILES:
            with open( fms.frag_sim_conf.file_path(fn), 'rb') as f:
                self.assertEqual( fms.__getattribute__(fn + '_hash'), hashlib.md5(f.read()).hexdigest())
        def add_mol():
             fms.frag_conf.create_frag_set(molecule = fms.molecule)
        self.assertRaises(Exception, add_mol, 'Unicity couple of molecule / frag_conf')

    def test_fs_gen_mgf(self):
        smiles = FragCommonTests.test_data['smiles']
        m = Molecule.load_from_smiles(smiles)
        fs = FragCommonTests.new_frag_sim()
        fm = fs.frag_molecule(m, adduct='M+')
        #'PEPMASS=45.0447245801'
        #'MOLECULAR_FORMULA=CH5N2+'
        exp_res = '\n'.join([\
                        'BEGIN IONS',\
                        'PEPMASS=45.0447245801',\
                        'CHARGE=1+',\
                        'FILENAME=METWORK_' + str(m.id), \
                        'SMILES=N=C[NH3+]',\
                        'INCHI=InChI=1S/CH4N2/c2-1-3/h1H,(H3,2,3)/p+1',\
                        'INCHIAUX=PNKUSGQVOMIXLU-UHFFFAOYSA-O',\
                        'SCANS=' + str(m.id), \
                        #'18.033826 0.397148', 
                        #'28.018175 0.993644', 
                        #'45.044725 98.609208', 
                        '18.0338255500 0.3971478062', 
                        '28.0181754800 0.9936442241', 
                        '45.0447245800 98.6092079700', 
                        'END IONS', 
                        ]) + u'\n'
        self.assertEqual(fm.gen_mgf(energy = 0, decimal = 10), exp_res)

