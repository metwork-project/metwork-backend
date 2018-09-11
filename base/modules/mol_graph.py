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
                'id': ion.ion_id,
                'parent_mass': ion.parent_mass,
            }
        } for ion in self.frag_sample.ions_list()]

        cosine_matrix = self.frag_sample.cosine_matrix
        matrix_size = len(cosine_matrix)

        edges = []

        for i in range(matrix_size):
            for j in range( matrix_size - i - 1 ) :
                cosine = cosine_matrix[i][j]
                if  cosine > 0:
                    edges.append({
                        'group': 'edges',
                        'data': {
                            'id': "{0} -- {1}".format(i, j),
                            'source': i,
                            'target': j,
                            'cosine': cosine
                        }
                    })

        return nodes + edges
