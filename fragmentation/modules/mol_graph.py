# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os, json
from libmetgem.network_generation import generate_network

class MolGraph:

    def __init__(self, frag_sample):
        self.frag_sample = frag_sample

    def gen_molecular_network(self):

        nodes = [{
            'group': 'nodes',
            'data': {
                'id': id,
                'name': ion.parent_mass,
                'parent_mass': ion.parent_mass,
            }
        } for id, ion in enumerate(self.frag_sample.ions_list())]

        cosine_matrix = self.frag_sample.cosine_matrix
        matrix_size = len(cosine_matrix)

        edges = [
            {
                'group': 'edges',
                'data': {
                    'id': "{0} -- {1}".format(i[0], i[1]),
                    'source': i[0],
                    'target': i[1],
                    'cosine': i[3]
                }
            }
        for i in generate_network(cosine_matrix, self.frag_sample.mzs(), 10, 0.65).tolist()
        if i[0] != i[1]]

        return nodes + edges
