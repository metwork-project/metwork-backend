#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

while [ "$1" != "" ]; do
    case $1 in
        --log) LOG=1;;
        --detached) DETACHED=1;;
        --test) TEST=1;;
        --stop) STOP=1;;
        *) echo "Option $1 not recognized"; exit 1;;
    esac
    shift
done

ROOT_DIR=$METWORK_BACKEND_PATH
ENV_DIR=$ROOT_DIR/envs

cd $ROOT_DIR

METWORK_VERSION=$(cat $ROOT_DIR/VERSION)


CONFIG="METWORK_CONFIG=$ENV_DIR/common.env,$ENV_DIR/local.env"
CMD="celery worker -A metwork_backend -Q "

if [ $LOG ]
then
    PID_PATH=$ROOT_DIR/metwork_worker_log.pid
    CMD="$CMD log.$METWORK_VERSION -l INFO --concurrency=1"
else
    PID_PATH=$ROOT_DIR/metwork_worker.pid
    CMD="$CMD web.$METWORK_VERSION,run.$METWORK_VERSION -l INFO"
fi

if [ "$TEST" ]
then
    CONFIG="$CONFIG,$ENV_DIR/test.env,$ENV_DIR/test-worker.env"
fi

if [ "$DETACHED" ]
then
    CMD="$CMD -D --pidfile=$PID_PATH"
fi

if test -f "$PID_PATH"; then
    kill $(cat $PID_PATH)
    rm $PID_PATH
fi

if [ -z $STOP ]
then
    export $CONFIG
    $CMD
fi

