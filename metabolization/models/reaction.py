# coding: utf8
from __future__ import unicode_literals

from django.db import models

import os.path
import hashlib
import subprocess
import json
import sys
import re
from django.conf import settings
from base.models import Tag
from base.modules import FileManagement, RDKit, ChemDoodle, ChemDoodleJSONError
from django.contrib.postgres.fields import JSONField
from django.db.models import Count


class Reaction(FileManagement, models.Model):
    class JSONAPIMeta:
        resource_name = "reactions"

    REACTANTS_MAX = 2

    name = models.CharField(max_length=128, default="new_reaction", unique=True)
    description = models.CharField(max_length=255, default="", null=True, blank=True)
    tags = models.ManyToManyField(Tag, related_name="reaction_tags", default=None)
    user = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, db_index=True
    )
    reactants_number = models.SmallIntegerField(default=0)
    smarts = models.CharField(
        max_length=1024, default=None, null=True, blank=True
    )  # smarts used by rdkit method
    status_code = models.PositiveSmallIntegerField(default=0, db_index=True)
    chemdoodle_json = JSONField(default=None, null=True, blank=True)
    chemdoodle_json_error = models.CharField(
        max_length=128, default=None, null=True, blank=True
    )

    public = True

    class status:
        INIT = 0
        EDIT = 10
        VALID = 20
        ACTIVE = 30
        OBSOLETE = 40
        ERROR = 90

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        if self.id is not None:
            prev_status = Reaction.objects.get(id=self.id).status_code
        else:
            prev_status = Reaction.status.EDIT
        if self.status_code == Reaction.status.INIT:
            self.status_code = Reaction.status.EDIT
        if (
            self.status_code == Reaction.status.EDIT
            and prev_status != Reaction.status.VALID
        ):
            try:
                cd = ChemDoodle()
                smarts = cd.json_to_react(self.chemdoodle_json)
                self.smarts = smarts
                react = RDKit.reaction_from_smarts(smarts)
                self.chemdoodle_json = cd.react_to_json(react)
                self.chemdoodle_json_error = None
                if self.ready():
                    self.status_code = Reaction.status.VALID
                else:
                    self.chemdoodle_json_error = "error with RDKit"
            except ChemDoodleJSONError as error:
                self.chemdoodle_json_error = error
            except:
                self.chemdoodle_json_error = "unexpected error"
        self.reactants_number = self.get_reactants_number()
        super(Reaction, self).save(*args, **kwargs)
        return self

    def tags_list(self):
        return [tag.name for tag in self.tags.all()]

    def load_smarts(self, smarts):
        try:
            cd = ChemDoodle()
            react = RDKit.reaction_from_smarts(smarts)
            self.smarts = smarts
            self.chemdoodle_json = cd.react_to_json(react)
            self.status_code = Reaction.status.VALID
        except:
            self.status_code = Reaction.status.EDIT
        self.save()
        return self

    def user_id(self):
        return self.user.id

    def user_name(self):
        return self.user.username

    @classmethod
    def activated(cls):
        return cls.objects.filter(status_code=cls.status.ACTIVE)

    @classmethod
    def max_delta(cls):
        max = 0
        for r in cls.activated():
            rd = r.mass_delta()
            if rd is not None and rd > max:
                max = rd
        return max

    @classmethod
    def create_from_smarts(cls, smarts, name, user, description=None):
        smarts = RDKit.reaction_to_smarts(RDKit.reaction_from_smarts(smarts))
        r = cls(name=name, user=user, description=description, smarts=smarts)
        r.save()
        return r

    def has_no_project(self):
        from metabolization.models import ReactionsConf

        return ReactionsConf.objects.filter(reactions__in=[self]).count() == 0

    def image_path(self):
        return "/".join([self.item_path(), "image.svg"])

    def ready(self):
        try:
            cd = ChemDoodle()
            self.chemdoodle_json = cd.react_to_json(
                RDKit.reaction_from_smarts(self.smarts)
            )
            rx = self.react_rdkit()
            return rx.Validate() == (0, 0)
        except:
            return False

    def react_rdkit(self):
        return self.react_rdkit_(self.smarts)

    def react_rdkit_(self, smarts):
        if smarts is not None:
            return RDKit.reaction_from_smarts(smarts)

    def get_reactants_number(self):
        smarts = self.smarts
        if smarts is not None:
            rx = self.react_rdkit_(smarts)
            return rx.GetNumReactantTemplates()
        else:
            return 0

    def run_reaction(self, reactants, method="rdkit"):
        if self.status_code < Reaction.status.ACTIVE:
            self.reactprocess_set.all().delete()
        from metabolization.models import ReactProcess

        rp = ReactProcess.objects.create()
        rp.reaction = self
        rp.reactants.set(reactants)
        rp.save()
        rp.run_reaction()
        return rp

    def mass_delta(self):
        from metabolization.models import ReactProcess

        rps = self.reactprocess_set.annotate(Count("products")).filter(
            products__count__gt=0
        )
        if len(rps) > 0:
            return rps.first().mass_delta()
        else:
            return None

    def update_name_tags(self):
        tag_pattern = re.compile(r"\[([\w\s]*)\]")
        name_pattern = re.compile(r"^([^\[\]]*\w)")

        new_name = re.findall(name_pattern, self.name)[0]
        tags = re.findall(tag_pattern, self.name)

        for tag_name in tags:
            tag_find = Tag.objects.filter(name=tag_name)
            if tag_find.count() == 1:
                tag = tag_find.first()
            else:
                tag = Tag.objects.create(name=tag_name)
            self.tags.add(tag)

        self.name = new_name
        self.save()
