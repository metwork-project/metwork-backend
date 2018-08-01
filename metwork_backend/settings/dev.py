from .defaults import *
import subprocess

DEBUG = True

# Rdkit warnings
RDLogger.logger().setLevel(RDLogger.CRITICAL)

SERVICES_HOST = ["db", "cache", "broker"]

try:
	SERVICES_IP = {
		key:subprocess.check_output(  ["getent", "hosts", key]).decode().split(' ')[0]
		for key in SERVICES_HOST }
except:
	SERVICES_IP = {
		key:"127.0.0.1"
		for key in SERVICES_HOST }

DATABASES = {
	'default': {
		'ENGINE': 'django.db.backends.postgresql_psycopg2',
		'NAME': os.environ['METWORK_DB_NAME'],
		'USER': 'metwork',
		'PASSWORD': os.environ['METWORK_DB_PASSWORD'], #'*metwork@db*',
		'HOST': SERVICES_IP['db'],#'localhost',
		'PORT': '5432',
	}
}

CACHES = {
	'default': {
		'BACKEND': 'django.core.cache.backends.memcached.MemcachedCache',
		'LOCATION': SERVICES_IP['cache'] + ':11211',
	}
}

CELERY_BROKER_URL = 'pyamqp://metwork:' + os.environ['METWORK_BROKER_PASSWORD'] + '@' + SERVICES_IP['broker'] +'/metwork'
