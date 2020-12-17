DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_DIR="$( cd $DIR/.. >/dev/null 2>&1 && pwd )"
ENV_DIR=$ROOT_DIR/envs


METWORK_VERSION=$(cat $ROOT_DIR/VERSION)
. $DIR/set-env.sh

BASE_CONFIG="METWORK_CONFIG=$ENV_DIR/common.env,$ENV_DIR/local.env"
BASE_CMD="celery worker -A metwork_backend -Q "

LOG_CMD="$BASE_CMD log.$METWORK_VERSION -l INFO --concurrency=1"
LOG_PID_PATH=$ROOT_DIR/metwork_worker_log.pid

PID_PATH=$ROOT_DIR/metwork_worker.pid
if [ "$1" == "log" ]
then
    PID_PATH=$LOG_PID_PATH
    BASE_CMD=$LOG_CMD
else
    BASE_CMD="$BASE_CMD web.$METWORK_VERSION,run.$METWORK_VERSION -l INFO"
fi

for PID_PATH_ in $PID_PATH $LOG_PID_PATH
do
    if test -f "$PID_PATH_"; then
        kill $(cat $PID_PATH_)
        rm $PID_PATH_
    fi
done

if [ "$1" == "test" ]
then
    export $BASE_CONFIG,$ENV_DIR/test.env,$ENV_DIR/test-worker.env
    $BASE_CMD -D --pidfile=$PID_PATH
    export $BASE_CONFIG,$ENV_DIR/test.env
elif [ "$1" == "test-log" ]
then
    export $BASE_CONFIG,$ENV_DIR/test.env,$ENV_DIR/test-worker.env
    $LOG_CMD -D --pidfile=$LOG_PID_PATH
    export $BASE_CONFIG,$ENV_DIR/test.env
elif [ "${@: -1}" == "detached" ]
then
    export $BASE_CONFIG
    $BASE_CMD -D --pidfile=$PID_PATH
elif [ "${@: -1}" != "stop" ]
then
    $BASE_CMD
fi
