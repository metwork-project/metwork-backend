# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TransactionTestCase
from django.contrib.auth import get_user_model
from base.models import Molecule  # , Project
from fragmentation.models import (
    FragAnnotation,
    FragAnnotationDB,
    FragAnnotationCompare,
    FragMolSample,
    FragSample,
)
from base.modules import BaseTestManagement


class FragAnnotationModelTests(BaseTestManagement):
    def test_import_file_annotation(self):
        user = get_user_model().objects.create()
        # project = Project.objects.create(user = user)
        data = [
            {
                "file_format": "default",
                "sample_file_path": "fragmentation/tests/files/test_import_annotation.mgf",
                "annotation_file_path": "fragmentation/tests/files/para_annotation.csv",
                "expected_res": [
                    (
                        1,
                        "CN=C1NC2=C(C)C3=C(C=CC2=N1)N=C(N3)N(C)C",
                        "DNP",
                        "DOL18-M M+H",
                    ),
                    (
                        2,
                        "[H]\C(CCNC(N)=N)=C1/N=C(O)N(\C([H])=C(/[H])C2=CC=C(O)C=C2)C1=O",
                        "DNP",
                        "PWJ39-P M+H",
                    ),
                    (
                        17,
                        "[H]\C(=C(\[H])C1=CC=C(O)C=C1)N1C(O)=NC(CCCNC(N)=N)C1=O",
                        "DNP",
                        "NGQ85-F M+H",
                    ),
                ],
            },
            {
                "file_format": "GNPS",
                "sample_file_path": "fragmentation/tests/files/test_import_annotation.mgf",
                "annotation_file_path": "fragmentation/tests/files/import_GNPS.tsv",
                "expected_res": [
                    (
                        1,
                        "CN=C1NC2=C(C)C3=C(C=CC2=N1)N=C(N3)N(C)C",
                        "GNPS : Commercial standard, Prasad",
                        "895096",
                    ),
                    (
                        36,
                        "C1=C(NC(=C1Br)Br)C(=O)NCC=CC2=CN=C(N2)N",
                        "GNPS : Other, CASMI",
                        "34649-22-4",
                    ),
                ],
            },
        ]

        for d in data:
            d["user"] = user
            self.eval_import_file_annotation(**d)

    def eval_import_file_annotation(
        self, file_format, sample_file_path, annotation_file_path, expected_res, user
    ):

        # sample_file_path = 'fragmentation/tests/files/test_import_annotation.mgf'
        with open(sample_file_path, "rb") as fss:
            fs = FragSample.import_sample(fss, user, "name", "file_name")
            fs.wait_import_done()

        # annotation_file_path = 'fragmentation/tests/files/para_annotation.csv'
        with open(annotation_file_path, "rb") as fa:
            # fs.import_annotation_file(fa)
            fs.import_annotation_file(fa, file_format)

        for er in expected_res:
            m = Molecule.find_from_smiles(er[1])
            self.assertTrue(m)
            fms_search = FragMolSample.objects.filter(frag_sample=fs, ion_id=er[0])
            self.assertEqual(fms_search.count(), 1)
            fms = fms_search.first()
            fa_search = FragAnnotationDB.objects.filter(
                molecule=m, frag_mol_sample=fms, db_source=er[2], db_id=er[3],
            )
            self.assertEqual(fa_search.count(), 1, fa_search.count())
