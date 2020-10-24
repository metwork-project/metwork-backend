#! bash
METWORK_VERSION=$(cat /opt/metwork-backend/VERSION)
echo "METWORK VERSION $METWORK_VERSION"
# ENV_PATH=/srv/metwork/conf/dev.env

# DJANGO_SETTINGS_MODULE=$(cat $ENV_PATH | grep DJANGO_SETTINGS_MODULE)
# export $DJANGO_SETTINGS_MODULE

# METWORK_CONFIG=$(cat $ENV_PATH | grep METWORK_CONFIG)
# export $METWORK_CONFIG

/opt/conda/envs/metwork/bin/celery worker -A metwork_backend -Q web.$METWORK_VERSION,run.$METWORK_VERSION -l INFO
