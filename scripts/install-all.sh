#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

# Pull submodules :
# git submodule update --init --recursive

SCRIPT_DIR=$METWORK_BACKEND_PATH/scripts

echo "######                       ########"
echo "######   Conda virtual env   ########"
echo "######                       ########"

if [ -z $METWORK_ENV ]
then
    METWORK_ENV=metwork
fi

if conda env list | grep -qw $METWORK_ENV
then
    echo "Conda $METWORK_ENV env already exist"
else
    echo "Creating $METWORK_ENV conda env"
    conda env create -n $METWORK_ENV -f conda-env.yml
fi

echo "######                       ########"
echo "######  Python dependencies  ########"
echo "######                       ########"

. $SCRIPT_DIR/install-python-depencies.sh

echo "######                       ########"
echo "######     Init Dataset      ########"
echo "######                       ########"

. $SCRIPT_DIR/init-dataset.sh

echo "######                       ########"
echo "######        CFM-ID         ########"
echo "######                       ########"

. $SCRIPT_DIR/install-cfm_id.sh
