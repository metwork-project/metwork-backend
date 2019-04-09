import os
import re
from django.core.exceptions import ImproperlyConfigured

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_app_version():
    with open(os.environ['METWORK_BACKEND_PATH'] + '/VERSION') as f:
        return f.read().strip()

APP_VERSION=get_app_version()

configs = {}

for path in re.split(',', os.environ.get('METWORK_CONFIG')):
    with open(path, 'r') as f:
    #configs = { **configs ,**json.loads(f.read())}
        configs = {
            **configs ,
            **{data[0]:data[1] for data in [
                l.replace('\n','').split('=') for l in f.readlines() if l !='\n' ] }
        }

def get_env(setting, configs=configs):
    try:
        val = configs[setting]
        if val == 'True':
            val = True
        elif val == 'False':
            val = False
        return val
    except KeyError:
        error_msg = "ImproperlyConfigured: Set {0} environment      variable".format(setting)
        raise ImproperlyConfigured(error_msg)
