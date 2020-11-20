SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. $SCRIPT_DIR/set-env.sh
. $SCRIPT_DIR/run-docker.sh
PID_PATH=$DIR/metwork_worker.pid

if test -f "$PID_PATH"; then
    echo "stop previous worker"
    kill $(cat $PID_PATH)
    rm $PID_PATH
fi

. $SCRIPT_DIR/run-worker.sh

./manage.py runserver

if test -f "$PID_PATH"; then
    echo "stop worker"
    kill $(cat $PID_PATH)
    rm $PID_PATH
fi

. $SCRIPT_DIR/run-docker.sh stop