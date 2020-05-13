# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import User, SampleAnnotationProject


@admin.register(User)
class UserAdmin(admin.ModelAdmin):
    readonly_fields = ("date_joined", "last_login", "samples_count", "projects_count")
    list_display = ("username", "samples_count", "projects_count")
    exclude = ("password", "first_name", "last_name")

    fieldsets = (
        (None, {"fields": ("email", "username", "organization")}),
        (
            "Activity",
            {
                "fields": (
                    "date_joined",
                    "last_login",
                    "samples_count",
                    "projects_count",
                ),
            },
        ),
        (
            "Authorization",
            {
                "fields": (
                    "is_active",
                    "is_staff",
                    "is_superuser",
                    "groups",
                    "user_permissions",
                ),
            },
        ),
    )

    def samples_count(self, instance):
        return instance.fragsample_set.count()

    samples_count.short_description = "Samples count"

    def projects_count(self, instance):
        return instance.project_set.count()

    projects_count.short_description = "Projects count"


@admin.register(SampleAnnotationProject)
class SampleAnnotationProjectAdmin(admin.ModelAdmin):
    readonly_fields = ("molecules_all_count", "frag_sample")
    fields = ("name", "user", "status_code", "frag_sample", "molecules_all_count")
    list_display = ("name", "user", "frag_sample", "status_code", "molecules_all_count")
    list_filter = ("status_code", "user")

    def molecules_all_count(self, instance):
        return instance.molecules_all_count()

    molecules_all_count.short_description = "Mols generated"
