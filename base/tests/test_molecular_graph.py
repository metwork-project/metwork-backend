# coding=utf-8
from __future__ import unicode_literals

from django.test import TransactionTestCase, Client
from django.db.utils import IntegrityError
from fragmentation.tests.test_frag_sample import FragSampleTestUtils
from base.tests.test_graph import GraphTestUtils
from base.models import MolecularGraph
from base.modules import BaseTestManagement


class MolecularGraphTestUtils(GraphTestUtils, FragSampleTestUtils):

    GRAPH_MODEL = MolecularGraph

    def create_graph(self, data=None, **kwargs):
        if not data:
            data = self.DEFAULT_DATA
        frag_sample = kwargs.pop("frag_sample", None)
        if not frag_sample:
            frag_sample = self.create_frag_sample()
        return frag_sample.molecular_network


class MolecularGraphTests(BaseTestManagement, MolecularGraphTestUtils):
    def test_error_if_create_without_frag_sample(self):

        with self.assertRaises(IntegrityError):
            self.GRAPH_MODEL.objects.create(data=self.DEFAULT_DATA)

    def test_create(self):
        frag_sample = self.create_frag_sample()

        graph = self.create_graph(frag_sample=frag_sample)

        assert graph.frag_sample == frag_sample
        assert frag_sample.molecular_network == graph

    def test_delete_with_frag_sample(self):

        graph = self.create_graph()
        assert self.GRAPH_MODEL.objects.get(id=graph.id)

        graph.frag_sample.delete()
        with self.assertRaises(self.GRAPH_MODEL.DoesNotExist):
            self.GRAPH_MODEL.objects.get(id=graph.id)

    def test_gen_data(self):

        graph = self.create_graph(data=None)

        graph.frag_sample.gen_cosine_matrix()

        graph.gen_data()

        assert graph.data != self.DEFAULT_DATA
