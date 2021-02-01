import os
from pathlib import Path
from .utils import get_env, PROJECT_ROOT

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.11/howto/static-files/

STATIC_URL = "/static/"
STATIC_ROOT = os.path.join(PROJECT_ROOT, "static")

DATA_FILES_PATH = str(
    Path(get_env("METWORK_DATA_PATH"), get_env("METWORK_DATA_FILES_FOLDER"))
)

if "METWORK_EDIT_FILES" in os.environ:
    EDIT_FILES = os.environ["METWORK_EDIT_FILES"] == "True"
else:
    EDIT_FILES = True

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql_psycopg2",
        "NAME": get_env("METWORK_DB_NAME"),
        "USER": "metwork",
        "PASSWORD": get_env("METWORK_DB_PASSWORD"),
        "HOST": get_env("METWORK_DB_HOST"),
        "PORT": "5432",
    }
}

# Cache
CACHES = {
    "default": {
        "BACKEND": "django.core.cache.backends.memcached.MemcachedCache",
        "LOCATION": get_env("METWORK_CACHE_HOST") + ":11211",
    }
}
