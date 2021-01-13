#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

while [ "$1" != "" ]; do
    case $1 in
        --no-worker) NO_WORKER=1;;
        --no-docker) NO_DOCKER=1;;
    esac
    shift
done

SCRIPT_DIR="$METWORK_BACKEND_PATH/scripts"

if [ -z $NO_DOCKER ]
then
    $SCRIPT_DIR/run-docker.sh --detached
fi

if [ -z $NO_WORKER ]
then
    $SCRIPT_DIR/run-worker.sh --detached
    $SCRIPT_DIR/run-worker.sh --log --detached
fi

$METWORK_BACKEND_PATH/manage.py runserver

if [ -z $NO_WORKER ]
then
    $SCRIPT_DIR/run-worker.sh --stop
    $SCRIPT_DIR/run-worker.sh --log --stop
fi

if [ -z $NO_DOCKER ]
then
    $SCRIPT_DIR/run-docker.sh --stop
fi