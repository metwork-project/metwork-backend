DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_DIR="$( cd $DIR/.. >/dev/null 2>&1 && pwd )"
ENV_DIR=$ROOT_DIR/envs


METWORK_VERSION=$(cat $ROOT_DIR/VERSION)
. $DIR/set-env.sh

BASE_CONFIG="METWORK_CONFIG=$ENV_DIR/common.env,$ENV_DIR/local.env"
BASE_CMD="celery worker -A metwork_backend -Q "

PID_PATH=$ROOT_DIR/metwork_worker.pid
if [ "$1" == "log" ]
then
    PID_PATH=$ROOT_DIR/metwork_worker_log.pid
    BASE_CMD="$BASE_CMD log.$METWORK_VERSION -l INFO --concurrency=1"
else
    BASE_CMD="$BASE_CMD web.$METWORK_VERSION,run.$METWORK_VERSION -l INFO"
fi

if test -f "$PID_PATH"; then
    kill $(cat $PID_PATH)
    rm $PID_PATH
fi


if [ "$1" == "test" ]
then
    export $BASE_CONFIG,$ENV_DIR/test.env,$ENV_DIR/test-worker.env
    $BASE_CMD -D --pidfile=$PID_PATH
    export $BASE_CONFIG,$ENV_DIR/test.env
elif [ "${@: -1}" == "detached" ]
then
    export $BASE_CONFIG
    $BASE_CMD -D --pidfile=$PID_PATH
elif [ "${@: -1}" != "stop" ]
then
    $BASE_CMD
fi
