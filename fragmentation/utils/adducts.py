# -*- coding: utf-8 -*-
import pandas as pd

class AdductManager:

    ADDUCTS_COLUMNS = ['label', 'mass', 'smiles']
    ADDUCTS_VALUES = {
        'positive': [
            ('M+', 0.0, ''),
            ('M+H', 1.007, 'H+'),
            ('M+NH4', 18.033825553, 'NH4+'),
            ('M+Na', 22.990, 'Na+'),
            ('M+K', 38.964, 'K+'),
        ],
        'negative': [
            ('M-', 0.0, ''),
            ('M-H', -1.007, ''),
            ('M+Cl', -34.969, 'Cl-'),
            ('M+HCOO', -44.998, '[O-]C=O'),
            ('M+CH3COO', -59.014, 'CC([O-])=O'),
        ]
    }
    MASS_TOLERANCE = 0.01

    def __init__(self, ion_charge='positive'):
        adducts = pd.DataFrame(self.ADDUCTS_VALUES[ion_charge], columns=self.ADDUCTS_COLUMNS)
        adducts.set_index('label', inplace=True)
        self.adducts = adducts
            

    def get_adduct(self, mol, frag_mol):
        """Get adduct of ion"""
        delta = frag_mol.parent_mass - mol.mass_exact()
        adducts = self.adducts
        try:
            return adducts.index[abs(adducts['mass']-delta)<=self.MASS_TOLERANCE][0]
        except IndexError:
            return None
