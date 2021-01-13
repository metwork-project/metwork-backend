SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_DIR="$( cd $SCRIPT_DIR/.. >/dev/null 2>&1 && pwd )"

git submodule update --init --recursive

. $SCRIPT_DIR/set-env.sh
conda install -c r -c rdkit -c conda-forge -c bioconda --yes --file $ROOT_DIR/conda-requirements.txt
pip install -r $ROOT_DIR/pip-requirements.txt