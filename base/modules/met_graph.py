# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.db.models import Count
import os, json
from django.db.models import Q
from django.db.models import Max
from fragmentation.utils import AdductManager

class MetGraph:

    def __init__(self, project):
        from fragmentation.models import FragAnnotationDB
        self.project = project
        self.adducts_mass = AdductManager().adducts.mass
        self.mols_in = [
            fa.molecule.id for \
                fa in FragAnnotationDB.objects.filter(
                    frag_mol_sample__frag_sample=project.frag_sample)]
        self.rps = [rp for rp in project.react_processes \
            .annotate(Count('products')) \
            .filter(products__count__gt=0) \
            .filter(Q(products__in=project.molecules_matching()) | Q(products__in=self.mols_in))]
        self.rps_dic = {val[1].id : val[0] for val in enumerate(self.rps)}
        self.mols = list( \
            set([m for rp in self.rps for m in rp.reactants.all()]) |
            set([m for rp in self.rps for m in self.products(rp)]))

        self.fm_ids = [fms.id for fms in self.project.frag_sample.fragmolsample_set.all()]
        self.fm_sims = {
            fac.frag_mol_compare.frag_mols.exclude(id=fac.frag_mol_sample.id).first() \
                for fac in self.project.fragannotationcompare_set.all()}
        self.fms_mols_mass = {
            fms.molecule.id: self.adducts_mass[fms.adduct] \
                for fms in self.fm_sims}
        self.mols_mass = {
            mol.id: mol.mass_exact() + self.get_adduct_mass(mol) \
                for mol in self.mols}
        self.mols_dic = { val[1].id : val[0] + len(self.rps) for val in enumerate(self.mols) }
        self.public_mols_id = {
            fad.molecule.id \
                for fad in FragAnnotationDB.objects.all() \
                if fad.is_public()}

    def get_adduct_mass(self, mol):
        from fragmentation.models import FragAnnotation
        # try:
        fas = FragAnnotation.objects.filter(molecule=mol, frag_mol_sample__in=self.fm_ids)
        if fas.count() > 0:
            if fas.first().adduct() is not None:
                return self.adducts_mass[fas.first().adduct()]
        return self.fms_mols_mass[mol.id]
        # except:
        #     return 0

    def products(self, rp):
        return rp.products.filter( Q(id__in=self.project.molecules_matching()) | Q(id__in=self.mols_in) ).distinct()

    def get_cosine(self, molecule):
        from fragmentation.models import FragAnnotationCompare
        query = \
            FragAnnotationCompare.objects\
                .filter(molecule = molecule, project = self.project)\
                .distinct()
        return ' | '.join([ \
            '{0} : {1}'.format(
            fac.frag_mol_sample.ion_id,
            fac.frag_mol_compare.cosine) \
            for fac in query.all() ])

    def get_best_cosine(self, molecule):
        from fragmentation.models import FragAnnotationCompare
        query = \
            FragAnnotationCompare.objects\
                .filter(molecule = molecule, project = self.project)\
                .aggregate(Max('frag_mol_compare__cosine'))
        return query['frag_mol_compare__cosine__max']

    def get_annotation_type(self, molecule):
        if molecule in self.project.molecules_init():
            return 'init'
        if molecule.id in self.public_mols_id:
            return 'public'
        return 'proposal'

    def metabolization_network(self):

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
                'reactJSON': rp.reaction.chemdoodle_json,
            }
        } for rp in self.rps]

        nodes += [{
            'group': 'nodes',
            'data': {
                'id': node_id('mol', m),
                'name': str(round(self.mols_mass[m.id], 3)),
                # 'name':  str(round( m.mass_exact(), 3 )),
                'parent_mass': str(round(self.mols_mass[m.id], 3)),
                'nodeType': 'molecule',
                'annotationType': self.get_annotation_type(m),
                'smiles': m.smiles(),
                'molJSON': m.chemdoodle_json,
                'cosine': self.get_cosine(m),
                'best_cosine': self.get_best_cosine(m),
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
        } for rp in self.rps for    m in rp.reactants.all() ]

        links += [{
        'group': 'edges',
            'data': {
                'id': "{0} -- {1}".format(node_id('react', rp), node_id('mol', m) ),
                'source': node_id('react', rp),
                'target': node_id('mol', m),
            }
        } for rp in self.rps for    m in self.products(rp)    ]

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
            for rp in self.rps for    m in rp.reactants.all() ]
        links += [{
            "id": "{0} -- {1}".format(self.rps_dic[rp.id], self.mols_dic[m.id]),
            "source": self.rps_dic[rp.id],
            "target": self.mols_dic[m.id],
            "interaction": "out",
            "reversible": "false"}
            for rp in self.rps for m in self.products(rp) ]
        res = '"nodes":' + json.dumps(nodes) +    ',\n"links":' + json.dumps(links)

        with open(self.project.item_path() + '/metwork_metexplore.json', 'w') as fw:
            with open('metabolization/modules/templates/metexplore_header', 'r') as fr:
                fw.writelines(fr.readlines())
            fw.writelines(res)
            with open('metabolization/modules/templates/metexplore_footer', 'r') as fr:
                    fw.writelines(fr.readlines())
