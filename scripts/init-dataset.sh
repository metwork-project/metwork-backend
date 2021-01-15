#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

SCRIPT_DIR="$METWORK_BACKEND_PATH/scripts"
. $SCRIPT_DIR/set-env.sh

echo "Create files directory"
. $METWORK_BACKEND_PATH/envs/local.env
mkdir $METWORK_DATA_FILES_PATH

$SCRIPT_DIR/run-docker.sh --detached
echo "Waiting for DB to initialize ..."
until PGPASSWORD=METWORK_DB_PASSWORD psql -q -d metwork -U metwork -h metwork_db -c ";" 2> /dev/null
do
    sleep 0.5s
done

echo "Migrate DB"
$METWORK_BACKEND_PATH/manage.py migrate

if test -f $METWORK_BACKEND_PATH/reactions_export.tsv
then
    echo "Create reactions"
    echo "from manage_reactions import import_reactions; import_reactions()" | ./manage.py shell
fi

$SCRIPT_DIR/run-docker.sh --stop
