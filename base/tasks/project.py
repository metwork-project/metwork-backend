from __future__ import absolute_import, unicode_literals
from celery import shared_task
from django.conf import settings
from base.models import Project
from base.modules.cache_management import model_from_cache


@shared_task
def finish_run(project_id):
    project = model_from_cache(Project, project_id)
    queue = settings.CELERY_WEB_QUEUE
    log_project.s(project_id, {"message": "run finished"}).apply_async(queue=queue)
    project.finish_run()


@shared_task
def log_project(project_id, data):
    project = model_from_cache(Project, project_id)
    project.log_to_file(data)
