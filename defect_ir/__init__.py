"""Defect Intermediate Representation (IDR) — unified 2D/3D defect graph layer.
See schema.py for the contract, validators.py for the checker, uncertainty.py for
candidate slip-system confidence, and adapters/ to lower an IDR to an engine input.
"""

from . import schema, validators, uncertainty
from .schema import empty_idr, edge_label, FORMAT_VERSION
from .validators import validate_idr, assert_valid
