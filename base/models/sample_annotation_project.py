# -*- coding: utf-8 -*-

# The aims of a SampleAnnotationProject is to create annotation
# of a FragMolSample of a FragSample with a Molecule generated with
# Reactions

from __future__ import unicode_literals
from pathlib import Path
import re
from django.db.models import Q
from django.conf import settings
from django.db import models
from base.models import Project, Molecule
from metabolization.models import Reaction, ReactionsConf, ReactProcess
from fragmentation.models import (
    FragSample,
    FragSimConf,
    FragAnnotationDB,
    FragCompareConf,
    FragAnnotationCompare,
)


class SampleAnnotationProject(Project):

    reactions_conf = models.ForeignKey(
        ReactionsConf, on_delete=models.PROTECT, default=None, null=True
    )
    frag_sample = models.ForeignKey(
        FragSample, on_delete=models.PROTECT, default=None, null=True
    )
    frag_annotations_init = models.ManyToManyField(FragAnnotationDB)
    frag_sim_conf = models.ForeignKey(
        FragSimConf, on_delete=models.PROTECT, default=None, null=True
    )
    frag_compare_conf = models.ForeignKey(
        FragCompareConf, on_delete=models.PROTECT, default=None, null=True
    )
    depth_total = models.IntegerField(default=0)
    depth_last_match = models.IntegerField(default=0)
    react_processes = models.ManyToManyField(ReactProcess)

    REACTIONS_LIMIT = settings.METWORK_CONF["MET"]["reactions_limit"]
    DEPTH_LIMIT = settings.METWORK_CONF["MET"]["depth_limit"]
    CUSTOM_FRAG_PARAMS_FILENAME = {
        "param": "param_output.log",
        "conf": "param_config.txt",
    }

    class JSONAPIMeta:
        resource_name = "projects"

    def save(self, *args, **kwargs):
        #  Manage conf of project when its new
        # Manage which information can be modified depends of project status_code
        # ===> To be used in Project ?
        had_pk = self.pk != None
        if not had_pk:
            super(SampleAnnotationProject, self).save(*args, **kwargs)
            self.load_default_conf()
        else:
            prev_status = SampleAnnotationProject.objects.get(id=self.id).status_code
            if max(self.status_code, prev_status) < Project.status.QUEUE:
                self.update_status()
                super(SampleAnnotationProject, self).save(*args, **kwargs)

            # While running only status_code can change.
            elif Project.status.RUNNING in (
                self.status_code,
                prev_status,
            ) or Project.status.QUEUE in (self.status_code, prev_status):
                kwargs["update_fields"] = [
                    "status_code",
                ]
                super(SampleAnnotationProject, self).save(update_fields=["status_code"])
            else:
                super(SampleAnnotationProject, self).save(
                    update_fields=["public", "name", "description"]
                )

        self.refresh_from_db()
        return self

    def clone_project(self):
        fields = [
            "description",
            "depth_total",
            "depth_last_match",
            "reactions_conf",
            "frag_sim_conf",
            "frag_compare_conf",
            "frag_sample",
        ]
        re_find = re.findall("(.*)(?:\sCOPY)[\s]?([\d]*)", self.name)
        if len(re_find) > 0:
            if re_find[0][1] == "":
                number = 1
            else:
                number = int(re_find[0][1]) + 1
            name = re_find[0][0] + " COPY " + str(number)
        else:
            name = self.name + " COPY"
        clone = SampleAnnotationProject(name=name, user=self.user)
        clone.save()
        for f in fields:
            clone.__setattr__(f, self.__getattribute__(f))
        clone.save()
        reaction_ids = [
            r.id
            for r in self.reactions_conf.reactions.exclude(
                status_code=Reaction.status.OBSOLETE
            )
        ]
        clone.change_reactions(reaction_ids)
        for fai in self.frag_annotations_init.all():
            clone.frag_annotations_init.add(fai)
        clone.save()
        return clone

    def reactions(self):
        if self.reactions_conf != None:
            return self.reactions_conf.reactions.all()

    def reactions_not_selected(self):
        return (
            Reaction.objects.all()
            .exclude(id__in=self.reaction_ids())
            .filter(status_code=Reaction.status.ACTIVE)
        )

    def reaction_ids(self):
        if self.reactions_conf != None:
            return [r.id for r in self.reactions_conf.reactions.all()]
        else:
            return []

    def frag_annotations_init_not_selected(self):
        selected_ids = [fs.id for fs in self.frag_annotations_init.all()]
        return FragAnnotationDB.objects.filter(
            frag_mol_sample__frag_sample=self.frag_sample
        ).exclude(id__in=selected_ids)

    def update_frag_sample(self, fs):
        # Assigned a FragSample to the project
        self.frag_annotations_init.clear()
        self.frag_sample = fs
        # By default, all annotation of FragSample are selected for the project
        conf_charge_params = FragSimConf.params_for_ion_charge(fs.ion_charge)
        self.update_conf_params("frag_sim_conf", **conf_charge_params)
        self.add_all_annotations()
        self.save()
        return self

    def add_all_annotations(self):
        for fa in FragAnnotationDB.objects.filter(
            frag_mol_sample__frag_sample=self.frag_sample
        ):
            self.frag_annotations_init.add(fa)

    def update_status(self):
        # Manage which status can be grant on the project depend on its state
        has_frag_sample = self.frag_sample != None
        self.molecules.clear()
        if has_frag_sample:
            for fm in self.frag_annotations_init.all():
                self.molecules.add(fm.molecule)
        has_molecules = self.molecules.count() > 0
        has_confs = (
            self.reactions_conf != None
            and self.frag_sim_conf != None
            and self.frag_compare_conf != None
        )
        if has_confs:
            has_reactions = self.reactions_conf.reactions.count() > 0
        if has_frag_sample and has_molecules and has_confs and has_reactions:
            self.status_code = 1
        else:
            self.status_code = 0
        return self

    def annotation_init_ids(self):
        if self.reactions_conf != None:
            return [fa.id for fa in self.frag_annotations_init.all()]

    def remove_reactions(self):
        self.change_reactions([])
        return self

    def select_reactions_by_tag(self):
        reaction_ids = [
            r.id
            for r in Reaction.activated().filter(tags__in=self.frag_sample.tags.all())
        ]
        self.change_reactions(reaction_ids)
        return self

    def add_all(self, dataLabel):
        if dataLabel == "reactions":
            reaction_ids = [r.id for r in Reaction.activated()]
            self.change_reactions(reaction_ids)
        if dataLabel == "annotations":
            self.add_all_annotations()
        return self

    def add_items(self, dataLabel, item_ids):
        if dataLabel == "reactions":
            reaction_ids = self.reaction_ids() + item_ids
            self.change_reactions(reaction_ids)
        if dataLabel == "annotations":
            for id in item_ids:
                annot = FragAnnotationDB.objects.get(id=id)
                self.frag_annotations_init.add(annot)
        return self

    def remove_all(self, dataLabel):
        if dataLabel == "reactions":
            self.remove_reactions()
        if dataLabel == "annotations":
            self.frag_annotations_init.clear()
        return self

    def remove_item(self, dataLabel, item_ids):
        item_id = item_ids[0]
        if dataLabel == "reactions":
            reaction = Reaction.objects.get(id=item_id)
            reaction_ids = self.reaction_ids()
            reaction_ids.remove(item_id)
            self.change_reactions(reaction_ids, reaction=reaction, to_remove=True)
        if dataLabel == "annotations":
            fa = FragAnnotationDB.objects.get(id=item_id)
            self.frag_annotations_init.remove(fa)
        return self

    def change_reactions(self, reaction_ids, reaction=None, to_remove=False):
        from django.db.models import Count

        # Manage which ReactionConf to applied to project
        reactions = Reaction.objects.filter(id__in=reaction_ids)
        # Look for existing conf that match the critteria
        rcs = (
            ReactionsConf.objects.annotate(Count("reactions"))
            .filter(
                reactions__count=len(reaction_ids),
                method_priority=self.reactions_conf.method_priority,
            )
            .distinct()
        )
        for r in reactions:
            rcs = rcs.filter(reactions__in=[r])
        prev_rc = self.reactions_conf
        # Does previous conf is not associated with other project ?
        prev_rc_uniq = self.reactions_conf.sampleannotationproject_set.count() == 1
        # Does new conf to applies already exist ?
        new_rc_exist = rcs.count() > 0
        if new_rc_exist:
            rc = rcs.first()
            self.reactions_conf = rc
            self.save()
            if prev_rc_uniq and prev_rc != rc:
                prev_rc.delete()
        elif reaction is not None and prev_rc_uniq:
            if to_remove:
                self.reactions_conf.reactions.remove(reaction)
            else:
                self.reactions_conf.reactions.add(reaction)
            self.reactions_conf.save()
        else:
            rc = ReactionsConf.objects.create(
                method_priority=self.reactions_conf.method_priority
            )
            for r in reactions:
                rc.reactions.add(r)
            rc.save()
            self.reactions_conf = rc
        self.save()

    def run(self):
        from base.tasks import start_run

        self.save()
        if self.status_code == Project.status.READY:
            self.status_code = Project.status.QUEUE
            self.save()
            start_run.delay(self.id)
            # If needed to start new project in priority
            # start.apply_async(args = [self.id], queue="priority.high")
        return self

    def molecules_(self, scope):
        fas_scope = {
            "init": self.frag_annotations_init.all(),
            "matching": FragAnnotationCompare.objects.filter(
                project=self, frag_mol_compare__match=True
            ),
        }
        if scope == "init_and_matching":
            return Molecule.objects.filter(
                Q(fragannotation__in=fas_scope["init"])
                | Q(fragannotation__in=fas_scope["matching"])
            ).distinct()
        else:
            return Molecule.objects.filter(
                fragannotation__in=fas_scope[scope]
            ).distinct()

    def molecules_matching(self):
        return self.molecules_("matching")

    def molecules_init(self):
        return self.molecules_("init")

    def molecules_init_and_matching(self):
        return self.molecules_("init_and_matching")

    def molecules_matching_count(self):
        return self.molecules_matching().count()

    def molecules_all_count(self):
        return self.molecules.count()

    def ms1_not_init(self):
        # return a list of id and a numpy array of m/z
        from fragmentation.models import FragMolSample, FragAnnotationDB
        import numpy as np

        fm_init = [fa.ion_id() for fa in self.frag_annotations_init.all()]
        fms = (
            FragMolSample.objects.filter(frag_sample=self.frag_sample)
            .exclude(ion_id__in=fm_init)
            .order_by("parent_mass")
        )
        return ([fm.id for fm in fms], np.array([float(fm.parent_mass) for fm in fms]))

    def finish_run(self, *args, **kwargs):
        from django.core.cache import cache

        super(SampleAnnotationProject, self).finish_run(*args, **kwargs)

        # Clear all cache concerning the project
        cache.delete("project_ms1_not_init_" + str(self.id))
        cache.delete("project_molecules_all_" + str(self.id))

        # Send email to user
        message = """
The run of the project {0} is finished.\n
Link to project : {1}/#/projects/{2}""".format(
            self.name, settings.FRONTEND_URL, self.id
        )
        self.user.email_user(subject="MetWork run finished", message=message)

    def gen_all_molecules(self):
        file_path = self.item_path() + "/metwork_all_molecules.csv"
        query_set = self.molecules.all()
        Molecule.gen_molecules(file_path, query_set)

    def gen_annotations(self):
        fas = FragAnnotationCompare.objects.filter(
            project=self, frag_mol_compare__match=True
        ).order_by("-frag_mol_compare__cosine")
        res = {}
        for fa in fas:
            ion_id = fa.frag_mol_sample.ion_id
            if ion_id in res:
                res[ion_id]["smiles"].append(fa.molecule.smiles())
            else:
                res[ion_id] = {
                    "smiles": [fa.molecule.smiles()],
                    "best_cosine": fa.frag_mol_compare.cosine,
                }
        res = "\n".join(
            [
                ",".join(
                    [
                        str(ion_id),
                        "|".join(res[ion_id]["smiles"]),
                        str(res[ion_id]["best_cosine"]),
                        "metwork",
                        str(ion_id),
                    ]
                )
                for ion_id in res
            ]
        )
        with open(self.item_path() + "/metwork_annotations.csv", "w") as fw:
            fw.writelines(
                ",".join(["ion_id", "smiles", "best cosine", "source", "sourceId"])
                + "\n"
            )
            fw.writelines(res)

    def gen_annotations_details(self):
        fas = FragAnnotationCompare.objects.filter(
            project=self, frag_mol_compare__match=True
        ).order_by("-frag_mol_compare__cosine")
        res = "\n".join(
            [
                ",".join(
                    [
                        str(fa.frag_mol_sample.ion_id),
                        fa.molecule.smiles(),
                        str(fa.frag_mol_compare.cosine),
                        "metwork",
                        str(fa.frag_mol_sample.ion_id),
                    ]
                )
                for fa in fas
            ]
        )
        with open(self.item_path() + "/metwork_annotations_details.csv", "w") as fw:
            fw.writelines(
                ",".join(["ion_id", "smiles", "cosine", "source", "sourceId"]) + "\n"
            )
            fw.writelines(res)

    def gen_metexplore(self):
        from base.modules import MetGraph

        mg = MetGraph(self)
        mg.gen_metexplore()

    def metabolization_network(self):
        from base.modules import MetGraph

        mg = MetGraph(self)
        return mg.metabolization_network()

    def load_custom_frag_param_files(self, file_type, data):
        self.save_custom_frag_param_files(file_type, data)
        file_path = self._get_custom_frag_param_path(file_type)
        self.update_conf_params(
            "frag_sim_conf", **{file_type + "_path": str(file_path)}
        )

    def save_custom_frag_param_files(self, file_type, data):
        file_path = self._get_custom_frag_param_path(file_type)
        file_path.write_text(data)

    def delete_custom_frag_param_files(self, file_type):
        file_path = self._get_custom_frag_param_path(file_type)
        if file_path.exists:
            file_path.unlink()

    def _get_custom_frag_param_path(self, file_type):
        return Path(self.item_path()) / self.CUSTOM_FRAG_PARAMS_FILENAME[file_type]
