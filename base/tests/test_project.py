# coding=utf-8
"""
Module description

:file path: /test_project.py
:authors: Yann Beauxis
:date: 2019-05-09
:python_version: 3.6
:license: BSD (3-clause)
"""
from __future__ import unicode_literals

import json
from pathlib import Path
from django.test import TransactionTestCase
from django.contrib.auth import get_user_model
from rest_framework.test import force_authenticate
from rest_framework.test import APIRequestFactory
from base.models import SampleAnnotationProject
from base.views import ProjectViewSet
from base.modules import BaseTestManagement


class ProjectRunModelTests(BaseTestManagement):
    def test_public(self):
        user_public = get_user_model().objects.create(
            username="pulbic", email="public@test.com"
        )
        user_private = get_user_model().objects.create(
            username="private", email="private@test.com"
        )
        project_public = SampleAnnotationProject.objects.create(
            name="public", user=user_private, public=True
        )
        project_private = SampleAnnotationProject.objects.create(
            name="private", user=user_private, public=False
        )

        TEST_PATTERN = [
            (project_private, user_private, 200),
            (project_public, user_private, 200),
            (project_private, user_public, 403),
            (project_public, user_public, 200),
        ]

        for project, user, status_code in TEST_PATTERN:
            resp = self.get_project(project, user)
            self.assertEqual(resp.status_code, status_code)

    def get_project(self, project, user):
        factory = APIRequestFactory()
        request = factory.get("/projects", format="json")
        force_authenticate(request, user=user)
        return ProjectViewSet.as_view({"get": "retrieve"})(request, pk=project.id)

    def test_log_to_file(self):
        user = get_user_model().objects.create(
            username="username", email="user@test.com"
        )
        project = SampleAnnotationProject.objects.create(name="name", user=user)
        datas = [
            {"task": "first", "result": "right"},
            {"task": "second", "result": "wrong"},
        ]
        for data in datas:
            project.log_to_file(data)
        dest_path = Path(project.item_path(), "log.json")
        exp_text = "\n".join([json.dumps(data) for data in datas]) + "\n"
        assert dest_path.read_text() == exp_text
