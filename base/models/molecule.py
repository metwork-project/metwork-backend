# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models
from base.modules import RDKit, ChemDoodle
from django_rdkit import models
from django.contrib.postgres.fields import JSONField

class Molecule(models.Model):

    name = models.CharField(
        max_length=255,
        default='')
    inchi_key = models.CharField(
        max_length=27,
        default=None,
        unique = True,
        db_index = True)
    mol_rdkit = models.MolField(
        default=None)
    smiles_with_limit = models.CharField(
        max_length=255,
        default='',
        db_index = True)
    chemdoodle_json = JSONField(
        default=None,
        null= True,
        blank=True)

    def __str__(self):
        return self.smiles()

    def __init__(self, *args, **kwargs):
        res = super().__init__(*args, **kwargs)
        if self.chemdoodle_json is None:
            self.chemdoodle_json = ChemDoodle().mol_to_json(self)
            self.save()
        return res

    class JSONAPIMeta:
        resource_name = "molecules"

    def smiles(self, kekulize = False):
        if (not kekulize) and (self.smiles_with_limit != ''):
            return self.smiles_with_limit
        else:
            return RDKit.mol_to_smiles(self.mol_rdkit, kekulize)

    def mass_exact(self):
        return RDKit.mass_exact(self.mol_rdkit)

    def mass_exact_isotopes(self):
    # return a list of masses including Brome isotopes
        mass = self.mass_exact()
        br_count = self.smiles().count('Br')
        return [ mass + 1.9979535*(inc) for inc in range(br_count + 1) ]

    def major_tautomers(self):
        import subprocess
        tautomer_major = set([])
        tautomer_all = subprocess.check_output([
                'cxcalc', 'tautomers',
                self.smiles(),
                '-f','smiles'])\
            .decode().split('\n')
        for sm in tautomer_all[:-1]:
            tautomer_major_ = subprocess.check_output([
                    'cxcalc', 'majortautomer',
                    sm,
                    '-f','smiles'])\
                .decode().split('\n')
            tautomer_major = tautomer_major | {Molecule.load_from_smiles(sm) for sm in tautomer_major_ }
        tautomer_major = tautomer_major - {False}
        return list(tautomer_major)

    @classmethod
    def create_from_smiles(cls, sm):
        try:
            m_rdkit = RDKit.mol_from_smiles(sm)
            if m_rdkit:
                m_inchi_key =  RDKit.mol_to_inchi_key(m_rdkit)
                m = Molecule(mol_rdkit = m_rdkit, inchi_key = m_inchi_key)
                sm_can = RDKit.mol_to_smiles(m_rdkit)
                if len(sm_can) < 255:
                    m.smiles_with_limit = sm_can
                m.save()
                return m
            else:
                raise Exception
        except:
            return False

    @classmethod
    def find_from_smiles(cls, sm):
    # Find a molecule from smiles by using inchi_key
        try:
            filter_res = Molecule.objects.filter(smiles_with_limit = sm)
            if filter_res.count() > 0:
                return filter_res.first()
            filter_res = Molecule.objects.filter(smiles_with_limit =\
                            RDKit.mol_to_smiles(RDKit.mol_from_smiles(sm)))
            if filter_res.count() > 0:
                return filter_res.first()
            filter_res = Molecule.objects.filter(inchi_key = \
                            RDKit.mol_to_inchi_key(RDKit.mol_from_smiles(sm)))
            if filter_res.count() > 0:
                return filter_res.first()
            return False
        except:
            return False

    @classmethod
    def load_from_smiles(cls, sm):
    # Return a Molecule instance if already exist (use find_from_smiles)
    # Or create it if correct smiles. Return False if error
        if sm != '':
            m = Molecule.find_from_smiles(sm)
            if not m : m = Molecule.create_from_smiles(sm)
            if not m:
                return False
            else:
                return m
        else:
            return False

    @classmethod
    def load_from_rdkit(cls, mol):
    # Return a Molecule instance corresponding to the rdkit mol in input
    # create it if not exist
        #from rdkit import Chem
        #Chem.Kekulize(mol)
        m_smiles = RDKit.mol_to_smiles(mol)
        return Molecule.load_from_smiles(m_smiles)

    @classmethod
    def gen_molecules(cls, file_path, queryset):
        with open(file_path, 'w') as fw:
            fw.writelines(','.join([
                'molecule_id',
                'inchi_key',
                'smiles',
            ]) + '\n')
            for m in queryset.order_by('id'):
                fw.writelines(','.join([
                    str(m.id),
                    m.inchi_key,
                    m.smiles(),
            ]) + '\n')
