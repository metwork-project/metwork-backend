#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

SCRIPT_DIR="$METWORK_BACKEND_PATH/scripts"
. $SCRIPT_DIR/set-env.sh
. $SCRIPT_DIR/run-docker.sh
. $SCRIPT_DIR/run-docker.sh stop

. $SCRIPT_DIR/run-docker.sh

./manage.py migrate

. $SCRIPT_DIR/run-docker.sh stop