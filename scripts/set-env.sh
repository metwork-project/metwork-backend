#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR="$( cd $DIR/.. >/dev/null 2>&1 && pwd )"
ENV_DIR=$DIR/envs

. $DIR/.env
export METWORK_BACKEND_PATH=$DIR
export METWORK_CONFIG=$DIR/.env,$ENV_DIR/common.env,$ENV_DIR/local.env

if [ -z $METWORK_ENV ]
then
    METWORK_ENV=metwork
fi

eval "$(conda shell.bash hook)"
conda activate $METWORK_ENV
