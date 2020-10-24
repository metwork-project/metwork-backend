# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import logging
from django.contrib.postgres.fields import JSONField
from polymorphic.models import PolymorphicModel

logger = logging.getLogger("django")


class GraphGenerator:

    DEFAULT_DATA = []

    def gen_graph(self):
        return self.DEFAULT_DATA


class Graph(PolymorphicModel):

    data = JSONField(blank=True, null=True)

    GRAPH_GENERATOR = GraphGenerator

    def get_data(self, force=False):
        if (not self.data) or force:
            self.gen_data()
        return self.data

    def gen_data(self):
        logger.debug("Generate graph data : Begin")
        self.data = self.GRAPH_GENERATOR(**self.generator_params()).gen_graph()
        self.save()
        logger.debug("Generate graph data : Done")

    def generator_params(self):
        return {}
