#! bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
METWORK_VERSION=$(cat $DIR/VERSION)
echo "METWORK VERSION $METWORK_VERSION"

export METWORK_BACKEND_PATH=$DIR
export METWORK_CONFIG=$DIR/common.env,$DIR/local.env
