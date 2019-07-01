# Celery settings
from .utils import APP_VERSION, get_env

#: Only add pickle to this list if your broker is secured
#: from unwanted access (see userguide/security.html)
CELERY_ACCEPT_CONTENT = ['json']
CELERY_RESULT_BACKEND = 'django-cache'
CELERY_TASK_SERIALIZER = 'json'
CELERY_TRACK_STARTED = True
CELERY_WEB_QUEUE = 'web.' + APP_VERSION
CELERY_RUN_QUEUE = 'run.' + APP_VERSION
CELERY_TASK_DEFAULT_QUEUE = CELERY_RUN_QUEUE
CELERY_QUEUES = {
    CELERY_WEB_QUEUE:{
        "exchange": CELERY_WEB_QUEUE,
        "routing_key": CELERY_WEB_QUEUE},
    CELERY_RUN_QUEUE:{
        "exchange": CELERY_RUN_QUEUE,
        "routing_key": CELERY_RUN_QUEUE}}

CELERY_BROKER_URL = \
    'pyamqp://metwork:' + get_env('METWORK_BROKER_PASSWORD') \
    + '@' + get_env('METWORK_BROKER_HOST') + '/metwork'
