# coding=utf-8
from __future__ import unicode_literals

from django.contrib.auth import get_user_model
from django.test import TransactionTestCase
from django.db.utils import IntegrityError
from base.models import Graph
from base.models.graph import GraphGenerator
from django.test import Client
from django.conf import settings

# from django.test.utils import setup_test_environment

# setup_test_environment()


class GraphTestUtils:

    GRAPH_MODEL = Graph
    DEFAULT_DATA = GraphGenerator.DEFAULT_DATA

    def create_graph(self, data, **kwargs):
        return self.GRAPH_MODEL.objects.create(data=data, **kwargs)


class GraphTests(TransactionTestCase, GraphTestUtils):
    def test_no_error_if_create_without_data(self):

        graph = self.create_graph(data=None)
        assert not graph.data

    def test_create_with_data(self):

        data = self.DEFAULT_DATA

        graph = self.create_graph(data)

        assert graph.data == self.DEFAULT_DATA

    def test_gen_data(self):

        graph = self.create_graph(data=None)

        graph.gen_data()

        assert graph.data == {
            "version": settings.APP_VERSION,
            "data": self.DEFAULT_DATA,
        }

    def test_get_data(self):

        for data_in, data_out in (
            (None, self.DEFAULT_DATA),
            ({"version": settings.APP_VERSION, "data": ["foo"],}, ["foo"]),
            ({"version": "0.0.1", "data": ["bar"],}, self.DEFAULT_DATA),
        ):

            graph = self.create_graph(data=data_in)

            assert graph.get_data() == data_out
