"""
WSGI config for metwork project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/howto/deployment/wsgi/
"""

import os
import sys
from django.core.wsgi import get_wsgi_application

settings_param_prod = "metwork_backend.settings.production"

path_prod = '/opt/metwork-backend'

if path_prod not in sys.path:
	sys.path.append(path_prod)

if not "DJANGO_SETTINGS_MODULE" in os.environ:
	os.chdir(path_prod)

os.environ["DJANGO_SETTINGS_MODULE"] = settings_param_prod

application = get_wsgi_application()

