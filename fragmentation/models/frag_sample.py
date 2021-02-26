# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
from decimal import *
import time
import re
import numpy as np
from celery import group, chain
from django.db import models, IntegrityError
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from libmetgem.cosine import compute_distance_matrix
from libmetgem.mgf import filter_data
from base.models import Molecule, Array2DModel, Tag
from base.modules import FileManagement
from fragmentation.utils import AdductManager, AnnotationStatus


class FragSample(FileManagement, models.Model, AdductManager):
    class JSONAPIMeta:
        resource_name = "fragsamples"

    class Meta:
        ordering = ("name",)

    user = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.PROTECT, default=None
    )
    name = models.CharField(max_length=128, default="")
    file_name = models.CharField(max_length=255, default="")
    tags = models.ManyToManyField(Tag, related_name="fragsample_tags", default=None)
    ion_charge = models.CharField(max_length=16, default="positive")
    description = models.CharField(max_length=255, default="", null=True, blank=True)
    ions_total = models.PositiveSmallIntegerField(default=0, db_index=True)
    cosine_matrix = models.OneToOneField(
        Array2DModel, related_name="cosine_matrix", on_delete=models.CASCADE, null=True,
    )
    status_code = models.PositiveIntegerField(default=0, db_index=True)

    conf = settings.METWORK_CONF["FRAG"]

    class status:
        INIT = 0
        READY = 1
        RUNNING = 2
        DONE = 3
        ERROR = 99

    def __str__(self):
        return self.name

    def is_public(self):
        for p in self.sampleannotationproject_set.all():
            if p.public:
                return True

    def obsolete(self):
        return True in (fms.obsolete() for fms in self.fragmolsample_set.all())

    def tags_list(self):
        return [tag.name for tag in self.tags.all()]

    def has_no_project(self):
        return self.sampleannotationproject_set.count() == 0

    def ions_count(self):
        return self.fragmolsample_set.count()

    def annotations_count(self):
        from fragmentation.models import FragAnnotationDB

        return FragAnnotationDB.objects.filter(
            frag_mol_sample__frag_sample=self
        ).count()

    @classmethod
    def import_sample(
        cls,
        file_object,
        user,
        name="",
        file_name="data.mgf",
        description="",
        ion_charge="positive",
        energy=1,
        task=False,
    ):

        from fragmentation import tasks

        pattern = re.compile(r"BEGIN IONS\n([\w\W\n]*?)END IONS")
        data = file_object.read().decode("utf-8")
        data = data.replace("\r\n", "\n").replace("\r", "\n")
        ions = re.findall(pattern, data)

        total_ions = len(ions)
        ions_limit = int(cls.conf["ions_limit"])
        if total_ions > ions_limit:
            raise IntegrityError(
                "{0} ions max authorized, {1} in the sample.".format(
                    ions_limit, total_ions
                )
            )

        fs = FragSample.objects.create(
            user=user,
            name=name,
            file_name=file_name,
            description=description,
            ion_charge=ion_charge,
            status_code=1,
            ions_total=total_ions,
        )
        fs.status_code = 2
        fs.save()
        fs.gen_item()
        with open(os.path.join(fs.item_path(), fs.file_name), "w") as fw:
            fw.writelines(data)

        if task:
            queue = settings.CELERY_WEB_QUEUE
            g = group(tasks.import_ion.s(fs.id, ion, energy) for ion in ions)
            s = chain(
                g,
                tasks.set_ions_count.s(fs.id),
                tasks.gen_cosine_matrix.s(),
                tasks.get_molecular_network.s(),
                tasks.finalize_import.s(),
            )
            s.apply_async(queue=queue)
            return fs
        else:
            return fs.import_sample_sync(ions, energy)

    def import_sample_sync(self, ions, energy):

        ions_count = 0

        for ion in ions:
            ion_count = self.import_ion(ion, energy)
            ions_count += ion_count
        self.gen_cosine_matrix()
        self.get_molecular_network()

        return self.finalise_import()

    def finalise_import(self):
        self.status_code = 3
        self.save()
        return self

    def import_ion(self, ion, energy):
        from fragmentation.models import (
            FragMolSample,
            FragMolAttribute,
            FragMolSpectrum,
        )

        params = re.findall(r"([^\n]*)=([^\n]*)", ion, re.U)
        pat = r"(\d+\.\d+)*E(\d+)"

        def conv_E(m):
            if m.group(1) is None:
                return m.group(0)
            else:
                return str(float(m.group(1)) * 10 ** int(m.group(2)))

        ion = re.sub(pat, conv_E, ion)
        peaks = re.findall(r"([\d]+\.+[\deE]+)[\t\s]([\d]+\.+[\deE]+)", ion, re.U)
        has_pepmass = "PEPMASS" in [v[0] for v in params]
        has_id = "SCANS" in [v[0] for v in params]
        has_peaks = len(peaks) > 1

        if has_pepmass and has_id and has_peaks:
            fsm = FragMolSample.objects.create(frag_sample=self)
            p = 1
            for param in params:
                if param[0] == "PEPMASS":
                    fsm.parent_mass = float(param[1])
                    # fsm.mass = Decimal(av[1])
                    fsm.save()
                elif param[0] == "SCANS":
                    fsm.ion_id = int(param[1])
                    fsm.save()
                else:
                    FragMolAttribute.objects.create(
                        frag_mol=fsm, title=param[0], value=param[1], position=p
                    )
                p += 1
            FragMolSpectrum.objects.create(
                frag_mol=fsm,
                spectrum=[[float(peak[0]), float(peak[1])] for peak in peaks],
                energy=energy,
            )
            return 1
        else:
            return 0

    def ions_list(self):
        return self.fragmolsample_set.all().order_by("ion_id").distinct()

    def mzs(self):
        return [fms.parent_mass for fms in self.ions_list()]

    def gen_cosine_matrix(self):
        query = self.ions_list()
        try:
            cosine_matrix = compute_distance_matrix(
                [fms.parent_mass for fms in query],
                [
                    filter_data(
                        np.array(fms.fragmolspectrum_set.get(energy=1).spectrum),
                        fms.parent_mass,
                        0.0,
                        0.0,
                        0.0,
                        0,
                    )
                    for fms in query
                ],
                0.002,
                5,
            )
            cosine_matrix = Array2DModel.objects.create(value=cosine_matrix.tolist())
            self.cosine_matrix = cosine_matrix
            self.save()
        except Exception as ex:
            pass
            # raise ex
        return self

    def get_molecular_network(self, task=False, force=False):
        from base.models import MolecularGraph

        try:
            self.molecular_network
        except FragSample.molecular_network.RelatedObjectDoesNotExist:
            MolecularGraph.objects.create(frag_sample=self)
        return self.molecular_network.get_data(force=force, task=task)

    def wait_import_done(self, timeout=360):
        begin = time.time()
        while self.status_code == FragSample.status.RUNNING:
            time.sleep(0.5)
            if (time.time() - begin) > timeout:
                print("\n#### close due to timeout #####\n")
                return self
            else:
                self.refresh_from_db()
        # time.sleep(20)
        return self

    def import_annotation_file(self, file_object, file_format="default"):
        from fragmentation.models import FragMolSample, FragAnnotationDB

        fls = [l.decode("utf-8") for l in file_object.readlines()]
        errors = {}
        col_titles = fls[0].split("\t")
        for i, fl in enumerate(fls[1:]):
            try:
                if file_format == "default":
                    ion_id, name, smiles, db_source, db_id = fl.split("\n")[0].split(
                        ","
                    )
                elif file_format == "GNPS":
                    data = fl.split("\t")
                    ion_id = data[col_titles.index("#Scan#")]
                    name = data[col_titles.index("Compound_Name")]
                    smiles = data[col_titles.index("Smiles")]

                    db_source = "GNPS : {0}, {1}".format(
                        data[col_titles.index("Compound_Source")],
                        data[col_titles.index("Data_Collector")],
                    )
                    db_id = data[col_titles.index("CAS_Number")]

                if int(ion_id) > 0:

                    molecule = Molecule.load_from_smiles(smiles)
                    fms = self.fragmolsample_set.get(ion_id=ion_id)
                    self.add_annotation(
                        frag_mol_sample=fms,
                        molecule=molecule,
                        name=name,
                        db_source=db_source,
                        db_id=db_id,
                        status_id=AnnotationStatus.VALIDATED,
                    )

            except Exception as err:
                # print(err)
                errors[i] = {"err": str(err), "smiles": smiles}
        return {"success": "Annotations successfully imported", "errors": errors}

    def add_annotation_from_smiles(
        self,
        ion_id,
        smiles,
        name="",
        db_source="",
        db_id="",
        status_id=AnnotationStatus.UNDEFINED,
    ):

        molecule = Molecule.load_from_smiles(smiles)
        fms = self.fragmolsample_set.get(ion_id=int(ion_id))

        self.add_annotation(
            frag_mol_sample=fms,
            molecule=molecule,
            name=name,
            db_source=db_source,
            db_id=db_id,
            status_id=status_id,
        )

    def add_annotation(
        self,
        frag_mol_sample,
        molecule,
        name="",
        db_source="",
        db_id="",
        status_id=AnnotationStatus.UNDEFINED,
    ):
        from fragmentation.models import FragAnnotationDB

        fms = frag_mol_sample
        am = AdductManager(ion_charge=self.ion_charge)
        adduct = am.get_adduct(molecule, fms)

        if adduct is not None:
            if fms.adduct != adduct:
                fms.adduct = adduct
                fms.save()

            if fms.adduct == adduct:
                FragAnnotationDB.objects.create(
                    frag_mol_sample=fms,
                    molecule=molecule,
                    name=name,
                    db_source=db_source,
                    db_id=db_id,
                    status_id=status_id,
                )

    def gen_mgf(self, energy=None, decimal=6):
        from fragmentation.models import FragMolSample

        res = "\n".join(
            [
                fm.gen_mgf(energy)
                for fm in FragMolSample.objects.filter(frag_sample=self).order_by(
                    "ion_id"
                )
            ]
        )
        return res

    def gen_mgf_file(self):
        self.gen_item()
        file_path = os.path.join(self.item_path(), self.file_name)
        if not os.path.exists(file_path):
            data = self.gen_mgf()
            with open(file_path, "w") as fw:
                fw.writelines(data)
