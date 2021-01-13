#!/bin/bash
. "$( dirname "$0" )/set-env.sh"

wget https://metwork.pharmacie.parisdescartes.fr/vendor/cfm_id.tar.gz -O $METWORK_BACKEND_PATH/cfm_id.tar.gz
tar -xf $METWORK_BACKEND_PATH/cfm_id.tar.gz -C $METWORK_BACKEND_PATH/../..
rm  $METWORK_BACKEND_PATH/cfm_id.tar.gz