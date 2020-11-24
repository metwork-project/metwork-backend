SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. $SCRIPT_DIR/set-env.sh
. $SCRIPT_DIR/run-docker.sh
. $SCRIPT_DIR/run-docker.sh stop

. $SCRIPT_DIR/run-docker.sh

./manage.py migrate

. $SCRIPT_DIR/run-docker.sh stop