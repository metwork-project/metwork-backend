#! bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DOCKER_DIR="$DIR/../docker"
CMD="docker-compose -f $DOCKER_DIR/docker-compose.yml"

if [ "$1" = "stop" ]
then
    $CMD stop
else
    $CMD up -d
fi