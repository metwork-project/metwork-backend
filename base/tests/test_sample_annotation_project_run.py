# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import time
from django.test import TransactionTestCase
from base.models import *
from metabolization.models import *
from fragmentation.models import *
from fragmentation.modules import FragSim
from django.contrib.auth import get_user_model
from metabolization.modules import ReactionTestManagement
from django.core.cache import cache
from django.test import tag


class SampleAnnotationProjectRunModelTests(ReactionTestManagement):
    @tag("integration")
    def test_batch(self):
        print("\n###### Begin project run tests ######")

        for p in self.params:
            print("\n###### {0} ######".format(p["name"]))
            self.eval_annotation_project(**p)

    def test_manual_stop(self):
        """Test if project is stopped manually before and of all run"""

        def run_test():
            self.eval_specific_annotation_project(name="depth_1", interrupt=0.5)

        self.assertRaises(Exception, run_test)

    def eval_specific_annotation_project(self, name="adducts", interrupt=None):
        for p in self.params:
            if p["name"] == name:
                self.eval_annotation_project(interrupt=interrupt, **p)

    def eval_annotation_project(
        self,
        name,
        anno_file,
        smiles,
        expected_anno,
        not_expected_smiles,
        reactions_name,
        react_process_count_expected=None,
        sample_file_name="test_annotation_project.mgf",
        interrupt=None,
    ):

        for m in [
            FragMolSim,
            SampleAnnotationProject,
            FragSimConf,
            Reaction,
            DefaultConf,
        ]:
            m.objects.all().delete()
        cache.clear()

        u = get_user_model().objects.create(email=name + "@test.com")

        sample_folder = "fragmentation/tests/files/test_annotation_project/"
        sample_file_path = sample_folder + sample_file_name
        anno_file_path = sample_folder + anno_file + ".csv"

        # Reaction.reactions_update()
        # reactions = [ Reaction.objects.get( name = rn) \
        #    for rn in reactions_name]

        smarts = {
            "methylation": "[N,O:1]>>[*:1]-[#6]",
            "diels_alder": "[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-1",
            "test_adducts": "[#8:3]-[#6]-[#6](-[#6]-[#8:1])-[#8:2]>>[#8:3]-[#6]-[#6](-[#6]-[#8:1]-[#6](=,:[#8])-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6])-[#8:2]",
        }

        reactions = [
            Reaction.create_from_smarts(smarts=smarts[rn], name=rn, user=u)
            for rn in reactions_name
        ]

        print(reactions)

        with open(sample_file_path, "rb") as fss:
            fs = FragSample.import_sample(fss, u)

        # Import annotation
        with open(anno_file_path, "rb") as f_annot:
            fs.import_annotation_file(f_annot)

        p = SampleAnnotationProject.objects.create(
            user=u,
        )
        for r in reactions:
            p.reactions.add(r)
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
        if interrupt is not None:
            print("begin wait")
            time.sleep(interrupt)
            print("end wait")
            p.finish_run()
            time.sleep(5)
        p.wait_run_end()
        self.assertEqual(p.status_code, Project.status.DONE)

        print(p.molecules_matching())

        print(".".join([m.smiles() for m in p.molecules.all()]) + "\n")

        for m in p.molecules.all():
            facs = FragAnnotationCompare.objects.filter(project=p, molecule=m)
            if facs.count() > 0:
                print(m.smiles(), facs.first().frag_mol_compare.cosine)

        # fsim = FragSim(fsc)
        # for m in p.molecules.all():
        #    print fsim.frag_molecule(m).gen_mgf()

        # Check if molecules expected added to project molecules
        # and annotation match

        for anno_exp in expected_anno:

            # Check molecule existence
            m = Molecule.find_from_smiles(anno_exp[0])
            self.assertTrue(m, name + " - " + anno_exp[0])
            self.assertIn(m, p.molecules.all())

            # Check annotation correct
            fms_search = FragMolSample.objects.filter(
                frag_sample=p.frag_sample, ion_id=anno_exp[1]
            )
            self.assertEqual(fms_search.count(), 1)
            # project = p,
            # frag_mol_sample = fms_search.first()).first().molecule.smiles()

            fac_search = FragAnnotationCompare.objects.filter(
                project=p,
                molecule=m,
                frag_mol_compare__match=True,
                frag_mol_sample=fms_search.first(),
            )
            self.assertNotEqual(
                fac_search.count(), 0, "fac_search.count() = " + str(fac_search.count())
            )

        # Check non expected molecules (depth to high) has not been generateds

        for sm_nexp in not_expected_smiles:
            m = Molecule.load_from_smiles(sm_nexp)
            self.assertFalse(m in p.molecules.all())

        # print '###### {0} ######'.format(name)
        # print ReactProcess.objects.count()
        # print '\n'.join([ '{0} with : {1}'.format(
        #                        rp.reaction,
        #                        ', '.join([ m.smiles() for m in rp.reactants.all() ])  ) \
        #                    for rp in ReactProcess.objects.all() ])

        # Check number of ReactProcess created optimize (not redondant)

        if react_process_count_expected != None:
            self.assertEqual(ReactProcess.objects.count(), react_process_count_expected)

        ReactProcess.objects.all().delete()

    params = [
        {
            "name": "depth_1",
            "anno_file": "anno_1",
            "smiles": ["OC=CCO"],
            "reactions_name": ["methylation"],
            "expected_anno": [("COC=CCO", 2)],
            "not_expected_smiles": ["COC=CCOC"],
        },
        {
            "name": "adducts",
            "sample_file_name": "test_adducts.mgf",
            "anno_file": "test_adducts_annot",
            "smiles": ["OCC(OC/C=C(C)/CC/C=C(C)/CC/C=C(C)/CC/C=C(C)/C)CO"],
            "reactions_name": ["test_adducts"],
            "expected_anno": [
                (
                    "OCC(OC/C=C(C)/CC/C=C(C)/CC/C=C(C)/CC/C=C(C)/C)COC(CCCCCCCCCCCCCCC)=O",
                    979,
                )
            ],
            "not_expected_smiles": [],
        },
        {
            "name": "bi_reactants_reaction",
            "anno_file": "anno_2",
            "smiles": ["OC=CCO", "C=Cc1c[nH]c(N)n1", "C=CCC"],
            "reactions_name": ["methylation", "diels_alder"],
            "expected_anno": [
                ("COC=CCO", 2),
                ("CCC1CCC=C2N=C(N)NC12", 10),
                ("CCC1CC=C2N=C(N)NC2C1", 9),
            ],
            "not_expected_smiles": [],
        },
    ]
