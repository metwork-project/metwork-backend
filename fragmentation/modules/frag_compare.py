# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from decimal import *
from fragmentation.models import FragMolCompare, FragMolPeak
from scipy import spatial
from django.core.exceptions import FieldError
from libmetgem.mgf import filter_data
from libmetgem.cosine import cosine_score
import numpy as np

class FragCompare:

    def __init__(self, frag_compare_conf):
        self.frag_compare_conf = frag_compare_conf

    def compare_frag_mols(self, frag_mols):
        if len(frag_mols) != 2:
            raise FieldError
        energies = [
            [fms[0] for fms in fm.fragmolspectrum_set.values_list('energy').order_by('energy').distinct()]
            for  fm in frag_mols ]
        best_cosine = 0
        best_fragmatch = 0
        best_energies = False
        for e0 in energies[0]:
            for e1 in energies[1]:
                params = [ (frag_mols[0], e0) , (frag_mols[1], e1) ]
                cosine= self.get_cosine(params)
                if cosine > best_cosine:
                    best_cosine = cosine
                    #best_fragmatch = num_frag_match
                    best_energies = (e0, e1)
        if best_energies:
            best_energies = "{0}:{1},{2}:{3}".format(frag_mols[0].id,best_energies[0],frag_mols[1].id,best_energies[1])
        fmc = FragMolCompare.objects.create(
            frag_compare_conf = self.frag_compare_conf,
            #frag_mols = frag_mols,
            cosine = best_cosine,
            #num_frag_match= best_fragmatch,
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
        parent_mass = [ \
            p[0].parent_mass \
            for p in params ]
        spectrum =  [ \
            p[0].fragmolspectrum_set.get(energy=p[1]).spectrum \
            for p in params]
#filter_data(data, mz_parent, min_intensity, parent_filter_tolerance, matched_peaks_window,min_matched_peaks_search)
        spectrum = [
            filter_data(np.array(spectrum[i]),parent_mass[i],0,0,0,1) \
            for i in range(2)]
#cosine_score(spectrum1_mz, spectrum1_data, spectrum2_mz, spectrum2_data, mz_tolerance, min_matched_peaks)
        return cosine_score(parent_mass[0], spectrum[0], parent_mass[1], spectrum[1], 10, 1)
