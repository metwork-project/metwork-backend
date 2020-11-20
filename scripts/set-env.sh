DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR="$DIR/.."
# DIR=$(pwd)
ENV_DIR=$DIR/envs

# METWORK_VERSION=$(cat $DIR/VERSION)
# echo "METWORK VERSION $METWORK_VERSION"

export METWORK_BACKEND_PATH=$DIR
export METWORK_CONFIG=$ENV_DIR/common.env,$ENV_DIR/local.env

if [ !METWORK_ENV ]
then
    METWORK_ENV=metwork
fi

source activate $METWORK_ENV
