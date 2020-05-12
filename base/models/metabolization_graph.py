# -*- coding: utf-8 -*-

from django.db import models

from base.models import Graph, SampleAnnotationProject
from base.modules import MetGraph


class MetabolizationGraph(Graph):

    project = models.OneToOneField(
        SampleAnnotationProject,
        on_delete=models.CASCADE,
        related_name="metabolization_graph",
    )

    GRAPH_GENERATOR = MetGraph

    def generator_params(self):
        return {"project": self.project}
