
from .defaults import *
DEBUG = False

# Rdkit warnings
RDLogger.logger().setLevel(RDLogger.CRITICAL)

DATABASES = {
	'default': {
		'ENGINE': 'django.db.backends.postgresql_psycopg2',
		'NAME': os.environ['METWORK_DB_NAME'],
		'USER': 'metwork',
		'PASSWORD': os.environ['METWORK_DB_PASSWORD'],
		'HOST': os.environ['METWORK_DB_HOST'],
		'PORT': '5432',
	}
}

# Cache
CACHES = {
	'default': {
		'BACKEND': 'django.core.cache.backends.memcached.MemcachedCache',
		'LOCATION': os.environ['METWORK_CACHE_HOST'] + ':11211',
	}
}

CELERY_BROKER_URL = 'pyamqp://metwork:' + os.environ['METWORK_BROKER_PASSWORD'] + '@' + os.environ['METWORK_BROKER_HOST'] + '/metwork'

