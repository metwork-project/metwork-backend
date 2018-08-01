# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from decimal import *
from fragmentation.models import FragMolCompare, FragMolPeak
from scipy import spatial
from django.core.exceptions import FieldError

class FragCompare:

    def __init__(self, frag_compare_conf):
        self.frag_compare_conf = frag_compare_conf

    def compare_frag_mols(self, frag_mols):
        if len(frag_mols) != 2:
            raise FieldError
        energies = [ 
            [fmp[0] for fmp in fm.fragmolpeak_set.values_list('energy').order_by('energy').distinct()]
            for  fm in frag_mols ]
        best_cosine = 0
        best_fragmatch = 0
        best_energies = False
        for e0 in energies[0]:
            for e1 in energies[1]:
                params = [ (frag_mols[0], e0) , (frag_mols[1], e1) ] 
                cosine, num_frag_match = self.get_cosine(params)
                if cosine > best_cosine:
                    best_cosine = cosine
                    best_fragmatch = num_frag_match
                    best_energies = (e0, e1)
        if best_energies:
            best_energies = "{0}:{1},{2}:{3}".format(frag_mols[0].id,best_energies[0],frag_mols[1].id,best_energies[1])
        fmc = FragMolCompare.objects.create(
            frag_compare_conf = self.frag_compare_conf,
            #frag_mols = frag_mols,
            cosine = best_cosine,
            num_frag_match= best_fragmatch,
            energies = best_energies)
        for fml in frag_mols:
            fmc.frag_mols.add(fml)
        self.evaluate_match(fmc)
        return fmc

    def evaluate_match(self, frag_mol_compare):
        frag_mol_compare.match = frag_mol_compare.cosine >= self.frag_compare_conf.cosine_threshold
        frag_mol_compare.save()
        return frag_mol_compare

    def get_cosine(self, params):
        _vectors = [\
                [\
                    [ float(fp.mz), float(fp.intensity) ] \
                for fp in FragMolPeak.objects.filter(frag_mol = p[0], energy=p[1]) ]\
            for p in params ]
    
        _dims = { vd[0] for v in _vectors for vd in v }
        _dims = list(_dims)
        _dims.sort()

        dims_dic = {}
        dims = []
        d_prev = 0
        for d in _dims:
            d_diff = abs(1 - d_prev / d) * 10**6
            if d_diff <= self.frag_compare_conf.ppm_tolerance:
            #if abs(d - d_prev) <= 0.02:
                dims_dic[d] = d_prev
                d_prev = 0
            else:
                dims_dic[d] = d
                d_prev = d
                dims.append(d)

        vectors_dic = { d : [0,0] for d in dims }
        inc = 0
        for fm in _vectors:
            for fp in fm:
                vectors_dic[ dims_dic[fp[0]] ][inc] = fp[1]
            inc += 1

        vectors = [[],[]]
        num_frag_match = 0
        for d in vectors_dic:
            for inc in range(2):
                intensity_raw = vectors_dic[d][inc]
                vectors[inc].append(intensity_raw**(1))
            if vectors_dic[d][0] * vectors_dic[d][1] > 0 :
                num_frag_match += 1
        cos = round(1.0 - spatial.distance.cosine(vectors[0], vectors[1]),3)
        cosine = Decimal(str(cos))
        return cosine, num_frag_match
