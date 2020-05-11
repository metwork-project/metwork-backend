# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from pathlib import Path
from django.contrib.auth import get_user_model
from base.models import Project, SampleAnnotationProject
from metabolization.models import Reaction, ReactionsConf
from metabolization.modules import ReactionTestManagement
from fragmentation.models import FragSample, FragSimConf


class SampleAnnotationProjectConfModelTests(ReactionTestManagement):

    user = None

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
            ("metabolization", "ReactionsConf", "reactions_conf"),
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
            set([r.id for r in p.reactions_conf.reactions.all()]),
            set([r.id for r in Reaction.objects.all()]),
        )

    def test_keep_custom_conf_while_save(self):
        p = self.create_project()
        rc = ReactionsConf.objects.create()
        p.reactions_conf = rc
        p.save()
        assert p.reactions_conf == rc

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

    def test_status_ready(self):
        initial_name = "initial name"
        # self.import_file(reaction_name = "methylation", user = u)
        self.create_reacts([("methylation", "[N,O:1]>>[*:1]-[#6]")])["methylation"]
        p = self.create_project(name=initial_name)
        self.assertEqual(p.status_code, Project.status.INIT)
        fs = self.create_frag_sample()
        p.update_frag_sample(fs)
        p.save()
        self.assertEqual(p.status_code, Project.status.INIT)
        fs.add_annotation(1, "CCC")
        p.update_frag_sample(fs)
        p.save()
        self.assertEqual(p.molecules.count(), 1)
        self.assertEqual(p.status_code, Project.status.INIT)
        rc = ReactionsConf.objects.create()
        self.assertEqual(rc.reactions.count(), 0)
        p.reactions_conf = rc
        p.save()
        self.assertEqual(p.status_code, 0)
        rc.reactions.add(Reaction.objects.first())
        p.save()
        self.assertEqual(p.status_code, Project.status.READY)
        p.run()
        self.assertTrue(p.status_code > Project.status.READY)
        other_name = "other"
        p.name = initial_name + " modified"
        p.save()
        self.assertEqual(p.name, initial_name)
        p.status_code = Project.status.DONE
        p.save()
        self.assertEqual(p.status_code, Project.status.DONE)
        p.status_code = Project.status.READY
        p.save()
        self.assertEqual(p.status_code, Project.status.DONE)

    def clone_project(self):
        # Reaction.reactions_update()
        initial_name = "origin name"
        clone_suffix = " COPY"
        p = self.create_project(name=initial_name)
        self.import_file(reaction_name="methylation", user=self.user)
        fs = self.create_frag_sample()
        p.update_frag_sample(fs)
        fs.add_annotation(1, "CCC")
        p.update_frag_sample(fs)
        rc = ReactionsConf.objects.create()
        p.reactions_conf = rc
        rc.reactions.add(Reaction.objects.first())
        p.save()
        pc = p.clone_project()
        self.assertNotEqual(pc, p)
        self.assertEqual(pc.name, p.name + clone_suffix)
        fields = [
            "user",
            "description",
            "depth_total",
            "depth_last_match",
            "reactions_conf",
            "frag_sim_conf",
            "frag_compare_conf",
            "frag_sample",
        ]
        self.assertEqual(pc.user, p.user)
        for f in fields:
            self.assertEqual(getattr(pc, f), getattr(p, f))

        def fais(project):
            return {fai.id for fai in project.frag_annotations_init.all()}

        self.assertEqual(fais(pc), fais(p))

    def test_save_custom_frag_param_files(self):
        from django.conf import settings

        p = self.create_project()
        CUSTOM_FILENAMES = {
            "param": "param_output_custom.log",
            "conf": "param_config_custom.txt",
        }
        sample_file_path = Path("fragmentation") / "tests" / "files" / "frag_sim_conf"
        for file_type in CUSTOM_FILENAMES:
            file_path = sample_file_path / CUSTOM_FILENAMES[file_type]
            data = file_path.read_text()
            target_file_name = SampleAnnotationProject.CUSTOM_FRAG_PARAMS_FILENAME[
                file_type
            ]
            target_file_path = Path(p.item_path()) / target_file_name
            if target_file_path.exists():
                target_file_path.unlink()

            p.save_custom_frag_param_files(file_type, data)
            assert target_file_path.exists
            assert target_file_path.read_text() == data

            target_file_path.unlink()
