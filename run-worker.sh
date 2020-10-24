#! bash
METWORK_VERSION=$(cat VERSION)
echo "METWORK VERSION $METWORK_VERSION"

BASE_CONFIG="METWORK_CONFIG=/home/yann/Dev/MetWork/metwork-backend/common.env,/home/yann/Dev/MetWork/metwork-backend/local.env"
BASE_CMD="celery worker -A metwork_backend -Q web.$METWORK_VERSION,run.$METWORK_VERSION -l INFO"

if [ "$1" == "test" ]
then
    export $BASE_CONFIG,/home/yann/Dev/MetWork/metwork-backend/test.env
    $BASE_CMD -D --pidfile=./metwork_worker.pid
else
    export $BASE_CONFIG
    $BASE_CMD
fi
#--concurrency=1 