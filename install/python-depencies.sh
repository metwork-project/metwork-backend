DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
conda install -c r -c rdkit -c conda-forge -c bioconda --yes --file $DIR/conda-requirements.txt
pip install -r $DIR/pip-requirements.txt