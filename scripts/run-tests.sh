SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. $SCRIPT_DIR/set-env.sh
# . $SCRIPT_DIR/run-docker.sh
PID_PATH=$DIR/metwork_worker.pid

# Remove previous test files
source $ENV_DIR/test.env 
rm -R $METWORK_DATA_FILES_PATH
mkdir $METWORK_DATA_FILES_PATH

if test -f "$PID_PATH"; then
    echo "stop previous worker"
    kill $(cat $PID_PATH)
    rm $PID_PATH
fi

. $SCRIPT_DIR/run-worker.sh test

echo "run unit tests"
./manage.py test --exclude-tag integration --failfast base.tests fragmentation.tests metabolization.tests


if test -f "$PID_PATH"; then
    echo "restart worker"
    kill $(cat $PID_PATH)
    rm $PID_PATH
fi

. $SCRIPT_DIR/run-worker.sh test

echo "run integration tests"
./manage.py test --tag integration --failfast base.tests
kill $(cat metwork_worker.pid)

# . $SCRIPT_DIR/run-docker.sh stop