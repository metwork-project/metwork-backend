git submodule update --init --recursive

. "$( dirname "$0" )/set-env.sh"

poetry install
