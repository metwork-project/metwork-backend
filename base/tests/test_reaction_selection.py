# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from django.contrib.auth import get_user_model
from base.models import Molecule, SampleAnnotationProject
from metabolization.models import Reaction
from fragmentation.models import FragSample

class ReactionSelectionTests(TransactionTestCase):

    def test_mass_selection(self):
        u = get_user_model().objects.create(email = 'user@test.com')
        methylation = Reaction.objects.create(
            user=u,
            name = 'methylation',
            smarts = '[N,O:1]>>[*:1]-[#6]' )
        methylation.run_reaction([Molecule.load_from_smiles('CCO')])
        diels_alder_modified = Reaction.objects.create(
            user=u,
            name = 'diels_alder_modified',
            smarts = '[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4]-[H].[#6:5]=,:[#6:6]>>[#6]-[#6:2]-1=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-[#6:1]-1' )
        diels_alder_modified.run_reaction([
            Molecule.load_from_smiles('C=Cc1c[nH]c(N)n1'),
            Molecule.load_from_smiles('C=C')])
        for r in [methylation, diels_alder_modified]:
            r.status_code = Reaction.status.ACTIVE
            r.save()

        sample_folder = 'fragmentation/tests/files/test_annotation_project/'
        sample_file_path = sample_folder + 'test_annotation_project.mgf'
        anno_file_path = sample_folder + 'anno_1.csv'

        with open(sample_file_path, 'rb') as fss:
            fs = FragSample.import_sample(fss, u, energy=2)

        # Import annotation
        with open(anno_file_path, 'rb') as f_annot:
            fs.import_annotation_file(f_annot)

        p = SampleAnnotationProject.objects.create(\
            user = u,
            depth_total = 3,
            depth_last_match = 0)
        p.save()
        p.update_frag_sample(fs)
        p.select_reaction_by_mass()

        self.assertIn(methylation, p.reactions())
        self.assertNotIn(diels_alder_modified, p.reactions())
