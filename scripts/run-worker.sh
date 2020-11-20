DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. $DIR/set-env.sh
ENV_DIR=$DIR/envs

BASE_CONFIG="METWORK_CONFIG=$ENV_DIR/common.env,$ENV_DIR/local.env"
BASE_CMD="celery worker -A metwork_backend -Q web.$METWORK_VERSION,run.$METWORK_VERSION -l INFO"

if [ "$1" == "test" ]
then
    export $BASE_CONFIG,$ENV_DIR/test.env,$ENV_DIR/test-worker.env
    $BASE_CMD -D --pidfile=./metwork_worker.pid
    export $BASE_CONFIG,$ENV_DIR/test.env
else
    export $BASE_CONFIG
    $BASE_CMD
fi
#--concurrency=1 