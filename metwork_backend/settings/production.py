
from .defaults import *
DEBUG = False

# Rdkit warnings
RDLogger.logger().setLevel(RDLogger.CRITICAL)

DATABASES = {
	'default': {
		'ENGINE': 'django.db.backends.postgresql_psycopg2',
		'NAME': get_env('METWORK_DB_NAME'),
		'USER': 'metwork',
		'PASSWORD': get_env('METWORK_DB_PASSWORD'),
		'HOST': get_env('METWORK_DB_HOST'),
		'PORT': '5432',
	}
}

# Cache
CACHES = {
	'default': {
		'BACKEND': 'django.core.cache.backends.memcached.MemcachedCache',
		'LOCATION': get_env('METWORK_CACHE_HOST') + ':11211',
	}
}

CELERY_BROKER_URL = 'pyamqp://metwork:' + get_env('METWORK_BROKER_PASSWORD') + '@' + get_env('METWORK_BROKER_HOST') + '/metwork'

