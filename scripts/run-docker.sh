#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

while [ "$1" != "" ]; do
    case $1 in
        --detached) DETACHED=1;;
        --stop) STOP=1;;
        --recreate) RECREATE=1;;
    esac
    shift
done

DOCKER_DIR="$METWORK_BACKEND_PATH/docker"
CMD="docker-compose -f $DOCKER_DIR/docker-compose.yml"

if [ "$STOP" ]
then
    $CMD stop
    exit 0
fi

CMD="$CMD up"

if [ "$RECREATE" ]
then
    CMD="$CMD --force-recreate"
fi


if [ "$DETACHED" ]
then
    CMD="$CMD -d"
fi

$CMD