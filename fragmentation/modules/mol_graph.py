# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os, json


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

        edges = []

        for i in range(matrix_size):
            for j in range( i + 1, matrix_size ) :
                cosine = cosine_matrix[i][j]
                if  cosine > 0.9:
                    edges.append({
                        'group': 'edges',
                        'data': {
                            'id': "{0} -- {1}".format(i, j),
                            'source': i,
                            'target': j,
                            'cosine': cosine
                        }
                    })

        print(len(edges))

        return nodes + edges
