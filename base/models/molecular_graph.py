# -*- coding: utf-8 -*-

from django.db import models

from base.models import Graph
from fragmentation.models import FragSample
from fragmentation.modules import MolGraph
from fragmentation.tasks import get_molecular_network


class MolecularGraph(Graph):

    frag_sample = models.OneToOneField(
        FragSample, on_delete=models.CASCADE, related_name="molecular_network"
    )

    GRAPH_GENERATOR = MolGraph
    TASK_GENERATOR = get_molecular_network

    def generator_params(self):
        if self.frag_sample.cosine_matrix is None:
            self.frag_sample.gen_cosine_matrix()
        return {"frag_sample": self.frag_sample}

    def task_args(self):
        return [self.frag_sample.id]
