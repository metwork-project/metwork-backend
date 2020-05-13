from __future__ import absolute_import, unicode_literals
from celery import shared_task
from base.models import Project
from base.modules.cache_management import model_from_cache


@shared_task
def finish_run(project_id):
    p = model_from_cache(Project, project_id)
    p.finish_run()
