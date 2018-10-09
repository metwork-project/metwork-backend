# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from base.models import *
from metabolization.models import *
from fragmentation.models import *
from fragmentation.modules import FragSim
from django.contrib.auth import get_user_model
from metabolization.modules import ReactionTestManagement
from django.core.cache import cache

class SampleAnnotationProjectRunModelTests(ReactionTestManagement):

    def test_batch(self):
        print ('\n###### Begin project run tests ######')
        for p in self.params:
            print ('\n###### {0} ######'.format(p['name']))
            self.eval_annotation_project(**p)

    def eval_specific_annotation_project(self, name = 'bromination'):
        for p in self.params:
            if p['name'] == name:
                self.eval_annotation_project(**p)

    def eval_annotation_project(self,
                name, anno_file,
                smiles, expected_anno,
                not_expected_smiles, reactions_name,
                depth_total, depth_last_match = 0,
                react_process_count_expected = None,
                sample_file_name = 'test_annotation_project.mgf'):

        for m in [FragMolSim, SampleAnnotationProject, FragSimConf, Reaction, DefaultConf]:
            m.objects.all().delete()
        cache.clear()

        u = get_user_model().objects.create(email = name + '@test.com')

        sample_folder = 'fragmentation/tests/files/test_annotation_project/'
        sample_file_path = sample_folder + sample_file_name
        anno_file_path = sample_folder + anno_file + '.csv'

        #Reaction.reactions_update()
        #reactions = [ Reaction.objects.get( name = rn) \
        #    for rn in reactions_name]

        smarts = {
            'methylation': '[N,O:1]>>[*:1]-[#6]',
            'diels_alder': '[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4]-[H].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-1',
        }

        reactions = [
            Reaction.create_from_smarts( smarts=smarts[rn], name=rn, user=u ) \
            for rn in reactions_name]

        print (reactions)

        rc = ReactionsConf.objects.create()
        for r in reactions:
            rc.reactions.add(r)


        with open(sample_file_path, 'rb') as fss:
            fs = FragSample.import_sample(fss, u, energy=2)

        # Import annotation
        with open(anno_file_path, 'rb') as f_annot:
            fs.import_annotation_file(f_annot)

        p = SampleAnnotationProject.objects.create(\
            user = u,
            #frag_sample = fs,
            depth_total = depth_total,
            depth_last_match = depth_last_match)
        p.reactions_conf = rc
        p.save()
        self.assertEqual(p.status_code, Project.status.INIT)

        p.update_frag_sample(fs)

        for sm in smiles:
            m = Molecule.find_from_smiles(sm)
            self.assertTrue(m)
            self.assertIn(m, p.molecules.all())

        # Check status and run

        p.update_status()
        self.assertEqual(p.status_code, Project.status.READY)
        p.run()
        self.assertTrue(p.status_code >= Project.status.QUEUE, p.status_code)
        p.wait_run_end()
        self.assertEqual(p.status_code, Project.status.DONE)

        print (p.molecules_matching())

        print ('.'.join([m.smiles() for m in p.molecules.all()]) + '\n')

        for m in  p.molecules.all():
            facs = FragAnnotationCompare.objects.filter(
                project = p,
                molecule = m)
            if facs.count() > 0:
              print (m.smiles(), facs.first().frag_mol_compare.cosine)

        #fsim = FragSim(fsc)
        #for m in p.molecules.all():
        #    print fsim.frag_molecule(m).gen_mgf()

        # Check if molecules expected added to project molecules
        # and annotation match

        for anno_exp in expected_anno:

            # Check molecule existence
            m = Molecule.find_from_smiles(anno_exp[0])
            self.assertTrue(m, name + ' - ' + anno_exp[0])
            self.assertIn(m, p.molecules.all())

            # Check annotation correct
            fms_search = FragMolSample.objects.filter(
                frag_sample = p.frag_sample,
                ion_id = anno_exp[1])
            self.assertEqual(fms_search.count(), 1)
                #project = p,
                #frag_mol_sample = fms_search.first()).first().molecule.smiles()

            fac_search = FragAnnotationCompare.objects.filter(
                project = p,
                molecule = m,
                frag_mol_compare__match=True,
                frag_mol_sample = fms_search.first())
            self.assertNotEqual(fac_search.count(), 0, 'fac_search.count() = ' + str(fac_search.count()))

        # Check non expected molecules (depth to high) has not been generateds

        for sm_nexp in not_expected_smiles:
            m = Molecule.load_from_smiles(sm_nexp)
            self.assertFalse(m in p.molecules.all())

        #print '###### {0} ######'.format(name)
        #print ReactProcess.objects.count()
        #print '\n'.join([ '{0} with : {1}'.format(
        #                        rp.reaction,
        #                        ', '.join([ m.smiles() for m in rp.reactants.all() ])  ) \
        #                    for rp in ReactProcess.objects.all() ])

        # Check number of ReactProcess created optimize (not redondant)

        if react_process_count_expected != None:
            self.assertEqual(ReactProcess.objects.count(), react_process_count_expected)

        ReactProcess.objects.all().delete()

    params = [ {
            'name' : 'depth_0',
            'anno_file' : 'anno_1',
            'smiles' : ['OC=CCO'],
            'reactions_name' : [ 'methylation' ],
            'expected_anno' : [],
            'not_expected_smiles' : ['COC=CCO'],
            'depth_total' : 0
        }, {
            'name' : 'depth_1',
            'anno_file' : 'anno_1',
            'smiles' : ['OC=CCO'],
            'reactions_name' : [ 'methylation' ],
            'expected_anno' : [('COC=CCO',2)],
            'not_expected_smiles' : ['COC=CCOC'],
            'depth_total' : 1
        }, {
            'name' : 'depth_2',
            'anno_file' : 'anno_1',
            'smiles' : ['OC=CCO'],
            'reactions_name' : [ 'methylation' ],
            'expected_anno' : [
                ('COC=CCO',2),
                ('COCC=CO', 3),
                ('COC=CCOC',15)],
            'not_expected_smiles' : [],
            'depth_total' : 2,
            'react_process_count_expected' : 3
        }, {
            'name' : 'depth_3',
            'anno_file' : 'anno_1',
            'smiles' : ['OC=CCO'],
            'reactions_name' : [ 'methylation' ],
            'expected_anno' : [
                ('COC=CCO',2),
                ('COCC=CO', 3),
                ('COC=CCOC',15)],
            'not_expected_smiles' : [],
            'depth_total' : 3,
            'react_process_count_expected' : 4
        }, {
            'name' : 'depth_2_gap_0',
            'anno_file' : 'anno_1',
            'smiles' : ['OC=CCO'],
            'reactions_name' : [ 'methylation' ],
            'expected_anno' : [],
            'not_expected_smiles' : ['COC=CCOC'],
            'depth_total' : 2,
            'depth_last_match' : 0,
            'sample_file_name' : 'test_annotation_project_GAP.mgf',
        }, {
            'name' : 'depth_2_gap_1',
            'anno_file' : 'anno_1',
            'smiles' : ['OC=CCO'],
            'reactions_name' : [ 'methylation' ],
            'expected_anno' : [('COC=CCOC',15)],
            'not_expected_smiles' : [],
            'depth_total' : 2,
            'depth_last_match' : 1,
            'sample_file_name' : 'test_annotation_project_GAP.mgf',
        }, {
            'name' : 'bi_reactants_reaction',
            'anno_file' : 'anno_2',
            'smiles' :['OC=CCO', 'C=Cc1c[nH]c(N)n1', 'C=CCC'],
            'reactions_name' : [ 'methylation', 'diels_alder' ],
            'expected_anno' : [
                ('COC=CCO',2),
                ('CCC1CCC=C2N=C(N)NC12', 10),
                ('CCC1CC=C2N=C(N)NC2C1', 9)],
            'not_expected_smiles' : [],
            'depth_total' : 1
        #}, {
        #    'name' : 'bromination',
        #    'anno_file' : 'para_annotation_single',
        #    'smiles' :[],
        #    'reactions_name' : [ 'bromination_of_phenols', 'bromination_of_phenols_isotope_81' ],
        #    'expected_anno' : [
        #        ('NC(=N)NCC\\C=C1\\N=C(O)N(\\C=C\\c2cc(Br)c(O)c([81Br])c2)C1=O',451),
        #        ],
        #    'not_expected_smiles' : [],
        #    'depth_total' : 10,
        #    'sample_file_name' : 'paraz 3 ions.mgf'
        } ]
