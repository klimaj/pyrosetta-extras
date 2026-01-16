"""Miscellaneous utilities for PyRosettaCluster unit tests."""

__author__ = "Jason C. Klima"


import platform


WEST_MIRROR = False
if WEST_MIRROR:
    ROSETTACOMMONS_CONDA_CHANNEL = "https://conda.rosettacommons.org"
    PYROSETTA_FIND_LINKS_PATH = "https://west.rosettacommons.org/pyrosetta/quarterly/release.cxx11thread.serialization"
else:
    ROSETTACOMMONS_CONDA_CHANNEL = "https://conda.graylab.jhu.edu"
    PYROSETTA_FIND_LINKS_PATH = "https://graylab.jhu.edu/download/PyRosetta4/archive/release-quarterly/release.cxx11thread.serialization"


def detect_platform():
    """Detect system platform string used by GitHub Actions."""
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == "linux":
        plat = "linux-64" if "64" in machine else "linux-32"
    elif system == "darwin":
        plat = "osx-arm64" if "arm" in machine else "osx-64"
    elif system == "windows":
        plat = "win-64"
    else:
        raise RuntimeError(f"Unsupported platform: {system} ({machine})")

    return plat
