DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
wget https://metwork.pharmacie.parisdescartes.fr/vendor/cfm_id.tar.gz -O $DIR/cfm_id.tar.gz
tar -xf $DIR/cfm_id.tar.gz -C $DIR/../..
rm  $DIR/cfm_id.tar.gz