# coding=utf-8
"""
Module description

:file path: /test_project.py
:authors: Yann Beauxis
:date: 2019-05-09
:python_version: 3.6
:license: BSD (3-clause)
"""
from __future__ import unicode_literals

from django.contrib.auth import get_user_model
from django.test import TransactionTestCase
from django.db.utils import IntegrityError
from base.models import GraphLayout
from django.test.utils import setup_test_environment
from django.test import Client


class GraphLayoutTests(TransactionTestCase):
    def create_user(self, user_email="user@test.com"):
        self.user = get_user_model().objects.create(email=user_email)

    def create_layout(
        self, graph_type=GraphLayout.types.FRAG_SAMPLE, data_id=1, layout={}
    ):

        return GraphLayout.objects.create(
            graph_type=graph_type, data_id=data_id, layout=layout
        )

    def test_create(self):

        with self.assertRaises(IntegrityError):
            layout = GraphLayout.objects.create()

        layout = {"key": "value"}

        graph_layout = self.create_layout(layout=layout)

        assert graph_layout.layout == layout

    def test_get_layout(self):
        params = {
            "data_id": 1,
            "graph_type": GraphLayout.types.FRAG_SAMPLE,
        }
        layout = {"key": "value"}
        self.create_layout(**params, layout=layout)
        client = Client()
        response = client.get("/graph-layouts", params)
        response = client.get("/tags")
        assert response.status_code == 200
