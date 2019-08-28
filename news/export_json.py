import json
from yaml import load, Loader
from django.conf import settings

news_conf = settings.METWORK_CONF['NEWS']

IMPORT_PATH = news_conf['yaml_path']
EXPORT_PATH = news_conf['json_path']
data = load(open(IMPORT_PATH, "r"), Loader=Loader)

with open(EXPORT_PATH, 'w') as fw:
    fw.write(json.dumps(data))