git submodule update --init --recursive

. "$( dirname "$0" )/set-env.sh"

conda install -c r -c rdkit -c conda-forge -c bioconda --yes --file $METWORK_BACKEND_PATH/conda-requirements.txt
pip install -r $METWORK_BACKEND_PATH/pip-requirements.txt