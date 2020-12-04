SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ENV_DIR="$( cd $DIR/../envs >/dev/null 2>&1 && pwd )"

. $SCRIPT_DIR/set-env.sh
. $SCRIPT_DIR/run-docker.sh
. $SCRIPT_DIR/run-worker.sh test

echo "run unit tests"
./manage.py test --exclude-tag integration --failfast base.tests fragmentation.tests metabolization.tests

. $SCRIPT_DIR/run-worker.sh stop
. $SCRIPT_DIR/run-worker.sh test

echo "run integration tests"
# Remove previous test files
source $ENV_DIR/test.env 
rm -R $METWORK_DATA_FILES_PATH
mkdir $METWORK_DATA_FILES_PATH

./manage.py test --tag integration --failfast base.tests


. $SCRIPT_DIR/run-worker.sh stop

# . $SCRIPT_DIR/run-docker.sh stop