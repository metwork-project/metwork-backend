# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import time
from pathlib import Path
from django.contrib.auth import get_user_model
from django.test import tag
from base.models import Project, SampleAnnotationProject
from metabolization.models import Reaction, ReactionsConf
from metabolization.modules import ReactionTestManagement
from fragmentation.models import FragSample, FragSimConf


class SampleAnnotationProjectConfModelTests(ReactionTestManagement):

    user = None
    CUSTOM_FRAG_FILENAMES = {
        "param": "param_output_custom.log",
        "conf": "param_config_custom.txt",
    }

    def create_user(self, user_email="user@test.com"):
        self.user = get_user_model().objects.create(email=user_email)

    def create_project(self, name="default_name"):
        if not self.user:
            self.create_user()
        return SampleAnnotationProject.objects.create(name=name, user=self.user,)

    def create_frag_sample(self, ion_charge="positive"):
        if not self.user:
            self.create_user()
        sample_file_path = "fragmentation/tests/files/example.mgf"
        with open(sample_file_path, "rb") as fss:
            fs = FragSample.import_sample(fss, self.user, "name", ion_charge=ion_charge)
            fs.wait_import_done()
        return fs

    def test_default_conf(self):
        p = self.create_project()
        for app_name, conf_class_name, p_conf_name in [
            ("fragmentation", "FragSimConf", "frag_sim_conf"),
            ("fragmentation", "FragCompareConf", "frag_compare_conf"),
        ]:
            # rc_filter = DefaultConf.objects.filter(
            #     project_class_name = project_class_name,
            #     app_name = app_name,
            #     conf_class_name = conf_class_name)
            # self.assertNotEqual( rc_filter.count(), 0)
            assert getattr(p, p_conf_name).id != 0
        self.assertEqual(
            set([r.id for r in p.reactions.all()]),
            set([r.id for r in Reaction.objects.all()]),
        )

    def test_change_params(self):
        p = self.create_project()
        CONF_NAME = "frag_sim_conf"
        DEFAULT_PATH = FragSimConf.PARAM_PATH_BY_CHARGE["positive"]
        PARAM_LABEL = "param_path"
        conf = getattr(p, CONF_NAME)
        assert getattr(conf, PARAM_LABEL) == DEFAULT_PATH
        NEW_PATH = FragSimConf.PARAM_PATH_BY_CHARGE["negative"]
        p.update_conf_params(CONF_NAME, **{PARAM_LABEL: NEW_PATH})
        conf = getattr(p, CONF_NAME)
        assert getattr(conf, PARAM_LABEL) == NEW_PATH

    def test_neg_params(self):
        for ion_charge in ("positive", "negative"):
            p = self.create_project()
            fs = self.create_frag_sample(ion_charge=ion_charge)
            p.update_frag_sample(fs)
            conf = getattr(p, "frag_sim_conf")
            assert (
                getattr(conf, "param_path")
                == FragSimConf.PARAM_PATH_BY_CHARGE[ion_charge]
            )
            assert (
                getattr(conf, "conf_path")
                == FragSimConf.CONF_PATH_BY_CHARGE[ion_charge]
            )

    @tag("integration")
    def test_status_ready(self):
        initial_name = "initial name"
        # self.import_file(reaction_name = "methylation", user = u)
        self.create_reacts([("methylation", "[N,O:1]>>[*:1]-[#6]")])
        p = self.create_project(name=initial_name)
        assert p.status_code == Project.status.INIT
        fs = self.create_frag_sample()
        p.update_frag_sample(fs)
        p.save()
        assert p.status_code == Project.status.INIT
        fs.add_annotation(1, "CCC")
        p.update_frag_sample(fs)
        p.save()
        assert p.molecules.count() == 1
        assert p.reactions.count() == 1
        assert p.status_code == Project.status.READY
        p.reactions.add(Reaction.objects.first())
        p.save()
        assert p.status_code == Project.status.READY
        p.run()
        self.assertTrue(p.status_code > Project.status.READY)
        other_name = "other"
        p.name = initial_name + " modified"
        p.save()
        assert p.name == initial_name
        p.status_code = Project.status.DONE
        p.save()
        assert p.status_code == Project.status.DONE
        p.status_code = Project.status.READY
        p.save()
        time.sleep(1)
        assert p.status_code == Project.status.DONE, p.status_code

    def test_clone_project(self):
        # Reaction.reactions_update()
        initial_name = "origin name"
        clone_suffix = " COPY"
        p = self.create_project(name=initial_name)
        self.create_reacts([("methylation", "[N,O:1]>>[*:1]-[#6]")])["methylation"]
        fs = self.create_frag_sample()
        p.update_frag_sample(fs)
        fs.add_annotation(1, "CCC")
        p.update_frag_sample(fs)
        p.reactions.add(Reaction.objects.first())
        p.save()
        assert len(p.reaction_ids()) == 1
        assert p.reactions.count() == 1
        pc = p.clone_project()
        assert pc != p
        assert pc.name == p.name + clone_suffix
        fields = [
            "user",
            "description",
            "depth_total",
            "depth_last_match",
            "frag_sim_conf",
            "frag_compare_conf",
            "frag_sample",
        ]
        assert pc.user == p.user
        for f in fields:
            assert getattr(pc, f) == getattr(p, f)

        def fais(project):
            return {fai.id for fai in project.frag_annotations_init.all()}

        assert fais(pc) == fais(p)

        assert p.reaction_ids() == pc.reaction_ids()

    def custom_frag_param_file_path(self, project, file_type):

        sample_file_path = Path("fragmentation") / "tests" / "files" / "frag_sim_conf"

        file_path = sample_file_path / self.CUSTOM_FRAG_FILENAMES[file_type]
        data = file_path.read_text()

        target_file_name = SampleAnnotationProject.CUSTOM_FRAG_PARAMS_FILENAME[
            file_type
        ]
        target_file_path = Path(project.item_path()) / target_file_name
        if target_file_path.exists():
            target_file_path.unlink()

        return target_file_path, data

    def test_save_custom_frag_param_files(self):

        project = self.create_project()

        for file_type in self.CUSTOM_FRAG_FILENAMES:
            target_file_path, data = self.custom_frag_param_file_path(
                project, file_type
            )

            project.save_custom_frag_param_files(file_type, data)
            assert target_file_path.exists()
            assert target_file_path.read_text() == data

            target_file_path.unlink()

    def test_delete_custom_frag_param_files(self):

        project = self.create_project()

        for file_type in self.CUSTOM_FRAG_FILENAMES:
            target_file_path, data = self.custom_frag_param_file_path(
                project, file_type
            )
            project.save_custom_frag_param_files(file_type, data)
            assert target_file_path.exists()
            project.delete_custom_frag_param_files(file_type)
            assert not target_file_path.exists()

    def test_load_custom_frag_param_files(self):

        project = self.create_project()

        for file_type in self.CUSTOM_FRAG_FILENAMES:
            target_file_path, data = self.custom_frag_param_file_path(
                project, file_type
            )

            project.load_custom_frag_param_files(file_type, data)
            assert target_file_path.exists()
            assert target_file_path.read_text() == data
            new_conf_path = getattr(project.frag_sim_conf, file_type + "_path")
            assert Path(new_conf_path) == target_file_path

    def test_reaction_conf_compatiblity(self):

        project = self.create_project()
        assert project.reactions_conf is None, project.reactions_conf

        self.create_reacts([("methylation", "[N,O:1]>>[*:1]-[#6]")])
        rc = ReactionsConf.objects.create()
        rc.reactions.add(Reaction.objects.first())
        project.reactions_conf = rc
        project.save()
        assert project.reactions_conf is not None
        reaction_ids = project.reaction_ids()
        assert len(reaction_ids) == 1

        clone = project.clone_project()
        assert clone.reactions_conf is None
        assert clone.reaction_ids() == reaction_ids

    def test_create_with_all_reactions(self):

        project = self.create_project("no reactions")
        assert project.reactions.count() == 0

        reacts = [
            ("methylation", "[N,O:1]>>[*:1]-[#6]"),
            (
                "diels_alder",
                "[#6:1]=,:[#6:2]-[#6:3]=,:[#6:4].[#6:5]=,:[#6:6]>>[#6:1]1-[#6:2]=,:[#6:3]-[#6:4]-[#6:6]-[#6:5]-1",
            ),
        ]
        self.create_reacts(reacts)

        project = self.create_project("with reactions")
        assert project.reactions.count() == len(reacts)

