from .base import *

DEBUG = True
FRONTEND_URL = 'http://localhost:4200'

ALLOWED_HOSTS += ['localhost']

EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
