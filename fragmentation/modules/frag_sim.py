# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import subprocess, re
from decimal import *
from fragmentation.models import FragSimConf, FragMol, FragMolSim, FragMolSpectrum
from django.db.models import Q
from fragmentation.utils import AdductManager


class FragSim:
    def __init__(self, frag_sim_conf):
        self.frag_sim_conf = frag_sim_conf

    def get_frag_mol(self, molecule, adduct="M+H"):
        # return FragMol.objects.instance_of(FragMolSim).filter(\
        fm_search = FragMolSim.objects.filter(
            frag_sim_conf=self.frag_sim_conf, molecule=molecule, adduct=adduct
        )
        if fm_search.count() > 0:
            return fm_search.first()
        else:
            return False

    def create_frag_mol_sim(self, molecule, adduct):
        fms = FragMolSim.objects.create(
            molecule=molecule,
            adduct=adduct,
            frag_sim_conf=self.frag_sim_conf,
            parent_mass=molecule.mass_exact(),
        )
        fms.update_hashes()
        return fms

    def frag_molecule(self, molecule, adduct="M+H", ion_charge="positive"):
        get_fm = self.get_frag_mol(molecule)
        if get_fm:
            return get_fm.wait_run_end()
        else:
            fms = self.create_frag_mol_sim(molecule, adduct)
            fms.status_code = FragMolSim.status.RUNNING
            fms.save()
            am = AdductManager(ion_charge)
            smiles = molecule.smiles()
            adduct_smiles = am.adducts.T[adduct].smiles
            if adduct_smiles != "":
                smiles = "{0}.[{1}]".format(smiles, adduct_smiles)
            run_out = subprocess.check_output(
                [
                    self.frag_sim_conf.CFM_ID_FOLDER + "cfm-predict",
                    smiles,
                    str(self.frag_sim_conf.threshold),
                    self.frag_sim_conf.file_path("param"),
                    self.frag_sim_conf.file_path("conf"),
                ]
            ).decode()
            _en, en = "", None
            spectrum = []
            first_en = True
            for sp in re.findall("energy(\d)\n([^energy]*)", run_out, re.U):
                FragMolSpectrum.objects.create(
                    frag_mol=fms,
                    energy=sp[0],
                    spectrum=[
                        [float(v) for v in peak.split(" ") if v != ""]
                        for peak in sp[1].split("\n")
                        if peak != ""
                    ],
                )

            fms.status_code = FragMolSim.status.DONE
            fms.save()
            return fms
