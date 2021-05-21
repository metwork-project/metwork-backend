# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import logging
import json
from django.contrib.postgres.fields import JSONField
from django.core.serializers.json import DjangoJSONEncoder
from polymorphic.models import PolymorphicModel
from django.conf import settings
from base.models import BaseModel

logger = logging.getLogger("django")


class GraphGenerator:

    DEFAULT_DATA = []

    def gen_graph(self):
        return self.DEFAULT_DATA


class Graph(PolymorphicModel, BaseModel):

    data = JSONField(blank=True, null=True)

    GRAPH_GENERATOR = GraphGenerator
    TASK_GENERATOR = None
    TASK_QUEUE = settings.CELERY_WEB_QUEUE

    def get_data(self, force=False, task=False):
        graph = getattr(self, "data", {})
        if (
            not isinstance(graph, dict)
            or graph.get("version", None) != settings.APP_VERSION
            or force
        ):
            return self.gen_data(task=task)
        return self.data["data"]

    def gen_data(self, task=False):
        if task and self.TASK_GENERATOR is not None:
            result = self.TASK_GENERATOR.s(*self.task_args()).apply_async(
                queue=self.TASK_QUEUE
            )
            result.get()
            self.refresh_from_db()
        else:
            logger.debug("Generate graph data : Begin")
            data = self.GRAPH_GENERATOR(**self.generator_params()).gen_graph()
            data = json.loads(json.dumps(data, cls=DjangoJSONEncoder))
            data = {"version": settings.APP_VERSION, "data": data}
            logger.debug("Generate graph data : Done")
            self.data = data
            self.save()
        return self.data["data"]

    def generator_params(self):
        return {}

    def task_args(self):
        return []
