#! bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
METWORK_VERSION=$(cat $DIR/VERSION)
echo "METWORK VERSION $METWORK_VERSION"

BASE_CONFIG="METWORK_CONFIG=$DIR/common.env,$DIR/local.env"
BASE_CMD="celery worker -A metwork_backend -Q web.$METWORK_VERSION,run.$METWORK_VERSION -l INFO"

if [ "$1" == "test" ]
then
    export $BASE_CONFIG,$DIR/test.env
    $BASE_CMD -D --pidfile=./metwork_worker.pid
else
    export $BASE_CONFIG
    $BASE_CMD
fi
#--concurrency=1 