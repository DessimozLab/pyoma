"""Linkout subpackage: generate crosslinks for external services pointing to OMA.

Submodules:
  ncbi    — NCBI LinkOut XML files and FTP upload
  uniprot — UniProtKB crosslink mapping
"""

from .ncbi import prepare_linkout_files, copy_to_ncbi

__all__ = [
    "prepare_linkout_files",
    "copy_to_ncbi",
]
