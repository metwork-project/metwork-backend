# coding=utf-8
from __future__ import unicode_literals
from django.test import TransactionTestCase, Client
from django.db.utils import IntegrityError
from django.contrib.auth import get_user_model
from fragmentation.tests.test_frag_sample import FragSampleTestUtils
from base.tests.test_graph import GraphTestUtils
from base.models import MetabolizationGraph, SampleAnnotationProject
from base.modules import BaseTestManagement


class MetabolizationGraphTestUtils(GraphTestUtils, FragSampleTestUtils):

    GRAPH_MODEL = MetabolizationGraph
    user = None

    def create_graph(self, data=None, **kwargs):
        if not data:
            data = self.DEFAULT_DATA
        project = kwargs.pop("project", None)
        if not project:
            frag_sample = self.create_frag_sample()
            project = self.create_project(frag_sample=frag_sample)
        return super().create_graph(data=data, project=project, **kwargs)

    def create_project(self, name="default_name", **kwargs):
        if not self.user:
            self.create_user()
        return SampleAnnotationProject.objects.create(
            name=name, user=self.user, **kwargs
        )

    def create_user(self, user_email="user@test.com"):
        self.user = get_user_model().objects.create(email=user_email)


class MetabolizationGraphTests(BaseTestManagement, MetabolizationGraphTestUtils):
    def test_error_if_create_without_frag_sample(self):

        with self.assertRaises(IntegrityError):
            self.GRAPH_MODEL.objects.create(data=self.DEFAULT_DATA)

    def test_create(self):
        project = self.create_project()

        graph = self.create_graph(project=project)

        assert graph.project == project
        assert project.metabolization_network == graph

    def test_delete_with_project(self):

        graph = self.create_graph()
        assert self.GRAPH_MODEL.objects.get(id=graph.id)

        graph.project.delete()
        with self.assertRaises(self.GRAPH_MODEL.DoesNotExist):
            self.GRAPH_MODEL.objects.get(id=graph.id)

    def test_gen_data(self):

        graph = self.create_graph(data=None)

        project = graph.project

        project.frag_sample.add_annotation(
            ion_id=1, smiles="CCC", db_source="None", db_id="1"
        )
        # project.run()
        # project.wait_run_end()

        graph.gen_data()

        # assert graph.data != self.DEFAULT_DATA
