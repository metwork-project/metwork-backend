SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_DIR="$( cd $SCRIPT_DIR/.. >/dev/null 2>&1 && pwd )"

. $SCRIPT_DIR/set-env.sh
. $SCRIPT_DIR/run-docker.sh
# . $SCRIPT_DIR/run-worker.sh detached

./manage.py runserver

. $SCRIPT_DIR/run-worker.sh stop
. $SCRIPT_DIR/run-docker.sh stop