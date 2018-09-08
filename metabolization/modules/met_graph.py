# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from base.models import *
from metabolization.models import *
from django.db.models import Count
import os, json
from django.db.models import Q

class MetGraph:

	def __init__(self, project):
		from fragmentation.models import FragAnnotationDB
		self.project = project
		self.mols_in = [fa.molecule.id  for fa in FragAnnotationDB.objects.filter(frag_mol_sample__frag_sample = project.frag_sample)]
		self.rps = [rp for rp in project.react_processes\
			.annotate(Count('products'))\
			.filter(products__count__gt = 0)\
			.filter( Q(products__in=project.molecules_matching()) | Q(products__in=self.mols_in) )]
		self.rps_dic = { val[1].id : val[0] for val in enumerate(self.rps) }
		self.mols = list( \
			set([ m for rp in self.rps for m in rp.reactants.all() ]) |
			set([ m for rp in self.rps for m in self.products(rp) ]) )
		self.mols_dic = { val[1].id : val[0] + len(self.rps) for val in enumerate(self.mols) }

	def products(self, rp):
		return rp.products.filter( Q(id__in=self.project.molecules_matching()) | Q(id__in=self.mols_in) ).distinct()

	def cytoscapejs_data(self):

		def node_id(node_type, element):
			if node_type == 'mol':
				dic_id = self.mols_dic[element.id]
			elif node_type == 'react':
				dic_id = self.rps_dic[element.id]
			return node_type + '_' + str(dic_id)

		nodes = [{
			'group': 'nodes',
			'data': {
				'id': node_id('react', rp),
				'name': rp.reaction.name,
				'nodeType': 'reaction',
				'reactionId': rp.reaction.id,
			}
		} for rp in self.rps]

		nodes += [{
			'group': 'nodes',
			'data': {
				'id': node_id('mol', m),
				'name': str(round( m.mass_exact(), 3 )) ,
				'nodeType': 'molecule',
				'annotation': 'init' if m in self.project.molecules_init() else 'proposal',
				'smiles': m.smiles(),
				'molFile': m.mol_file(),
			}
		} for m in self.mols]

		links = [{
			'group': 'edges',
			'data': {
				'id': "{0} -- {1}".format( node_id('mol', m), node_id('react', rp) ),
				'source': node_id('mol', m),
				'target': node_id('react', rp),
				#'source': 'mol' + self.mols_dic[m.id],
				#'target': 'react_' + self.rps_dic[rp.id],
			}
		} for rp in self.rps for  m in rp.reactants.all() ]

		links += [{
			'group': 'edges',
			'data': {
				'id': "{0} -- {1}".format(node_id('react', rp), node_id('mol', m) ),
				'source': node_id('react', rp),
				'target': node_id('mol', m),
			}
		} for rp in self.rps for  m in self.products(rp)  ]

		return nodes + links
		#return json.dumps(nodes + links )

	def gen_metexplore(self):
		res = ''

		nodes = [{'name': rp.reaction.name, "biologicalType":"reaction"} for rp in self.rps]
		nodes += [{'name': str(round(m.mass_exact(),3)), "biologicalType":"metabolite", "labelVisible":False} for m in self.mols]
		links = [{
			'id': "{0} -- {1}".format(self.mols_dic[m.id], self.rps_dic[rp.id]),
			'source': self.mols_dic[m.id],
			'target': self.rps_dic[rp.id],
			"interaction":"in",
			"reversible":"false"}
			for rp in self.rps for  m in rp.reactants.all() ]
		links += [{
				"id": "{0} -- {1}".format(self.rps_dic[rp.id], self.mols_dic[m.id]),
				"source": self.rps_dic[rp.id],
				"target": self.mols_dic[m.id],
				"interaction": "out",
				"reversible": "false"}
				for rp in self.rps for m in self.products(rp) ]
		res = '"nodes":' + json.dumps(nodes) +  ',\n"links":' + json.dumps(links)


		with open(self.project.item_path() + '/metexplore.json', 'w') as fw:
			with open('metabolization/modules/templates/metexplore_header', 'r') as fr:
				fw.writelines(fr.readlines())
			fw.writelines(res)
			with open('metabolization/modules/templates/metexplore_footer', 'r') as fr:
				fw.writelines(fr.readlines())
