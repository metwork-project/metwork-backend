import os
from .base import *

DEBUG = True
FRONTEND_URL = 'http://localhost:4200'

ALLOWED_HOSTS += ['localhost', 'testserver']

EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'


LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
        },
    },
    'loggers': {
        'django': {
            'handlers': ['console'],
            'level': os.getenv('DJANGO_LOG_LEVEL', 'INFO'),
        },
    },
}