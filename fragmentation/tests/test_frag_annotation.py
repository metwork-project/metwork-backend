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

from fragmentation.utils import LEVEL_STATUS_MAPPING

from base.modules import BaseTestManagement


class FragAnnotationModelTests(BaseTestManagement):
    TEST_DATA = {
        "level": {
            "file_format": "level",
            "sample_file_path": "fragmentation/tests/files/test_import_annotation.mgf",
            "annotation_file_path": "fragmentation/tests/files/para_annotation_level.csv",
            "expected_res": [
                (
                    "1",
                    "1",
                    "CN=C1N=C2C=CC3=C(NC(N(C)C)=N3)C(C)=C2N1",
                    "DNP",
                    "DOL18-M M+H",
                ),
                (
                    "2",
                    "2",
                    "N=C(N)NCC/C=C1/N=C(O)N(/C=C/c2ccc(O)cc2)C1=O",
                    "DNP",
                    "PWJ39-P M+H",
                ),
                (
                    "17",
                    "17",
                    "N=C(N)NCCCC1N=C(O)N(/C=C/c2ccc(O)cc2)C1=O",
                    "DNP",
                    "NGQ85-F M+H",
                ),
            ],
            "level": (30, 20, 10),
        },
        "GNPS": {
            "file_format": "GNPS",
            "sample_file_path": "fragmentation/tests/files/test_import_annotation.mgf",
            "annotation_file_path": "fragmentation/tests/files/import_GNPS.tsv",
            "expected_res": [
                (
                    "1",
                    "5-HYDROXY-L-TRYPTOPHAN",
                    "CN=C1NC2=C(C)C3=C(C=CC2=N1)N=C(N3)N(C)C",
                    "GNPS : Commercial standard, Prasad",
                    "895096",
                ),
                (
                    "36",
                    "Oroidin - CASMI2016 Category 1 - Challenge 2",
                    "C1=C(NC(=C1Br)Br)C(=O)NCC=CC2=CN=C(N2)N",
                    "GNPS : Other, CASMI",
                    "34649-22-4",
                ),
            ],
        },
    }

    def test_import_file_annotation(self):
        user = get_user_model().objects.create()

        data = self.TEST_DATA.copy()
        data["default"] = data["level"]
        data["default"].pop("level")
        data["default"]["file_format"] = "default"
        data["default"][
            "annotation_file_path"
        ] = "fragmentation/tests/files/para_annotation.csv"

        for d in data.values():
            d["user"] = user
            self.eval_import_file_annotation(**d)

    def test_gen_annotations_file(self):
        data = self.TEST_DATA["level"].copy()
        level = data.pop("level")
        expected_res = data.pop("expected_res")
        expected_res = [
            tuple((*res, LEVEL_STATUS_MAPPING[level[idx]]))
            for idx, res in enumerate(expected_res)
        ]
        user = get_user_model().objects.create()
        frag_sample = self.import_annotation_file(**data, user=user)
        data = frag_sample.gen_annotations()
        assert data == expected_res
        file_path = frag_sample.gen_annotations_file()
        expected_res = (
            "\n".join(
                [
                    ",".join([str(val) for val in row])
                    for idx, row in enumerate(expected_res)
                ]
            )
            + "\n"
        )
        assert file_path.read_text() == expected_res

    def import_annotation_file(
        self, file_format, sample_file_path, annotation_file_path, user,
    ):

        with open(sample_file_path, "rb") as fss:
            frag_sample = FragSample.import_sample(fss, user, "name", "file_name")
            frag_sample.wait_import_done()

        with open(annotation_file_path, "rb") as fa:
            frag_sample.import_annotation_file(fa, file_format)

        return frag_sample

    def eval_import_file_annotation(
        self,
        file_format,
        sample_file_path,
        annotation_file_path,
        expected_res,
        user,
        level=None,
    ):

        frag_sample = self.import_annotation_file(
            file_format, sample_file_path, annotation_file_path, user
        )

        for idx, er in enumerate(expected_res):
            m = Molecule.find_from_smiles(er[2])
            self.assertTrue(m)
            fms_search = FragMolSample.objects.filter(
                frag_sample=frag_sample, ion_id=int(er[0])
            )
            self.assertEqual(fms_search.count(), 1)
            fms = fms_search.first()
            fa_search = FragAnnotationDB.objects.filter(
                molecule=m,
                frag_mol_sample=fms,
                name=er[1],
                db_source=er[3],
                db_id=er[4],
            )
            self.assertEqual(fa_search.count(), 1, fa_search.count())
            if level:
                assert fa_search.first().status_id == level[idx]

