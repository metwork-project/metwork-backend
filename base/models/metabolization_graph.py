# -*- coding: utf-8 -*-

from django.db import models

from base.models import Graph, SampleAnnotationProject
from base.modules import MetGraph
from base.tasks import get_metabolization_network


class MetabolizationGraph(Graph):

    project = models.OneToOneField(
        SampleAnnotationProject,
        on_delete=models.CASCADE,
        related_name="metabolization_network",
    )

    GRAPH_GENERATOR = MetGraph
    TASK_GENERATOR = get_metabolization_network

    def generator_params(self):
        return {"project": self.project}

    def task_args(self):
        return [self.project.id]
