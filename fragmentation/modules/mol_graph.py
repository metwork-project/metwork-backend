# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os, json
from libmetgem.network_generation import generate_network
from fragmentation.models import FragAnnotationDB, FragAnnotationCompare
from fragmentation.utils import AnnotationStatus


class MolGraph:
    def __init__(self, frag_sample):
        self.frag_sample = frag_sample

    def annotation_status(self, ion):
        query = FragAnnotationDB.objects.filter(
            frag_mol_sample__frag_sample_id=self.frag_sample.id
        ).filter(frag_mol_sample__id=ion.id)
        for status in ("VALIDATED", "PUTATIVE", "UNRATED"):
            if query.filter(status_id=AnnotationStatus[status]).count() > 0:
                return AnnotationStatus[status].value
        return AnnotationStatus.UNDEFINED

    def best_annotation(self, ion):
        best_annot = ion.best_annotation()
        mol = best_annot[0]
        if mol is not None:
            smiles = mol.smiles()
        else:
            smiles = None

        return {"smiles": smiles, "cosine": best_annot[1]}

    def mol_json(self, ion):
        best_annot = ion.best_annotation()[0]
        if best_annot is not None:
            return best_annot.chemdoodle_json

    def gen_graph(self):

        nodes = [
            {
                "group": "nodes",
                "data": {
                    "id": id,
                    "name": ion.parent_mass,
                    "parent_mass": ion.parent_mass,
                    "nodeType": "ion",
                    "annotationStatusId": self.annotation_status(ion),
                    "bestAnnotation": self.best_annotation(ion),
                    "info": ion.gen_info().replace("\n", "<br/>"),
                    "molJSON": self.mol_json(ion),
                },
            }
            for id, ion in enumerate(self.frag_sample.ions_list())
        ]

        cosine_matrix = self.frag_sample.cosine_matrix.value
        matrix_size = len(cosine_matrix)

        edges = [
            {
                "group": "edges",
                "data": {
                    "id": "{0} -- {1}".format(i[0], i[1]),
                    "source": i[0],
                    "target": i[1],
                    "cosine": i[3],
                },
            }
            for i in generate_network(
                cosine_matrix, self.frag_sample.mzs(), 10, 0.65
            ).tolist()
            if i[0] != i[1]
        ]

        return nodes + edges
