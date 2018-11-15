# Chemistry
from rdkit import RDLogger
from .utils import get_env

# RDKit loggin level
RDLogger.logger().setLevel(RDLogger.CRITICAL)

CFM_ID_PATH=get_env("METWORK_CFM_ID_PATH")

PROTON_MASS = 1.007
