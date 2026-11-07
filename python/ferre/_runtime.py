"""Locate and run the native tools included in an installed FERRE wheel."""

import os
from pathlib import Path
import sys


def get_executable(name="ferre"):
    """Return the absolute path to an installed FERRE tool.

    Parameters
    ----------
    name : str, optional
        Either ferre or ascii2bin.

    Returns
    -------
    str
        Path to the executable bundled with this package.

    Raises
    ------
    ValueError
        The requested tool name is unsupported.
    FileNotFoundError
        This package has no compiled executable. Install a platform wheel
        or build a wheel from source with GNU Fortran.
    """
    if name not in ("ferre", "ascii2bin"):
        raise ValueError(f"Unknown FERRE tool: {name}")
    path = Path(__file__).resolve().parent / "bin" / name
    if not path.is_file():
        raise FileNotFoundError(
            f"Compiled {name} is missing. Install a platform wheel or build "
            "from source with GNU Fortran."
        )
    return str(path)


def main():
    """Forward command-line arguments and the working directory to FERRE."""
    os.execv(get_executable(), ["ferre", *sys.argv[1:]])


def ascii2bin_main():
    """Run the bundled ASCII-to-binary library converter."""
    os.execv(get_executable("ascii2bin"), ["ascii2bin", *sys.argv[1:]])
