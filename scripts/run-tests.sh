#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

TEST_LIST=(
    "--exclude-tag integration --failfast base.tests fragmentation.tests metabolization.tests"
    "--tag integration --failfast base.tests"
)

while [ "$1" != "" ]; do
    case $1 in
        --unit) TEST_LIST=("${TEST_LIST[0]}");;
        --integration) TEST_LIST=("${TEST_LIST[1]}");;
        --no-worker) NO_WORKER=1;;
        --no-docker) NO_DOCKER=1;;
        --recreate-docker) RECREATE_DOCKER="--recreate";;
        *) TEST_LIST=("$1")
        # *) echo "Option $1 not recognized"; exit 1;;
    esac
    shift
done

cd $METWORK_BACKEND_PATH

SCRIPT_DIR="$METWORK_BACKEND_PATH/scripts"
ENV_DIR="$METWORK_BACKEND_PATH/envs"
METWORK_DATA_FILES_PATH=$METWORK_DATA_PATH/files_test

if [ -z $NO_DOCKER ]
then
    $SCRIPT_DIR/run-docker.sh --detached $RECREATE_DOCKER
fi

source $ENV_DIR/test.env 

for params in "${TEST_LIST[@]}"; do

    rm -R $METWORK_DATA_FILES_PATH
    mkdir $METWORK_DATA_FILES_PATH

    if [ -z $NO_WORKER ]
    then
        $SCRIPT_DIR/run-worker.sh --test --detached
        $SCRIPT_DIR/run-worker.sh --test --log --detached
    fi

   ./manage.py test $params

    if [ -z $NO_WORKER ]
    then
    $SCRIPT_DIR/run-worker.sh --test --stop
    $SCRIPT_DIR/run-worker.sh --test --log --stop
    fi

done

if [ $METWORK_TEST_SOUND ]
then
    paplay $METWORK_TEST_SOUND
fi