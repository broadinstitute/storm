"""Build TR rule sidecar tables (Parquet) for Hail joins."""

from .builder import SidecarResult, build_tr_sidecar
from .schema import RULESET_VERSION, SCHEMA_VERSION

__all__ = [
    "RULESET_VERSION",
    "SCHEMA_VERSION",
    "SidecarResult",
    "build_tr_sidecar",
]
