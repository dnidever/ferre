"""Build FERRE executables inside a platform-specific Python wheel."""

import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py
from wheel.bdist_wheel import bdist_wheel


class BuildFerre(build_py):
    """Compile in a temporary directory so source checkouts remain clean."""

    def run(self):
        super().run()
        source = Path(__file__).resolve().parent / "src"
        makefile = (source / "makefile").read_text().replace("\\\n", " ")
        modules = re.search(r"^SOURCE_MODULES\s*=\s*(.*)$", makefile, re.M)[1].split()
        routines = re.search(r"^SOURCE\s*=\s*(.*)$", makefile, re.M)[1].split()
        compiler = shlex.split(os.environ.get("FC", "gfortran"))
        flags = shlex.split(os.environ.get(
            "FFLAGS", "-O2 -fopenmp -fallow-argument-mismatch"
        ))
        link_flags = shlex.split(os.environ.get("LDFLAGS", ""))
        destination = Path(self.build_lib).resolve() / "ferre" / "bin"
        destination.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="ferre-build-") as temporary:
            work = Path(temporary)
            for path in source.glob("*.f90"):
                shutil.copy2(path, work)
            for name, sources in (
                ("ferre", modules + routines),
                ("ascii2bin", ["share.f90", "timer.f90", "checkstat.f90",
                              "quit.f90", "ascii2bin.f90"]),
            ):
                subprocess.run(compiler + flags + sources + link_flags + ["-o", name],
                               cwd=work, check=True)
                target = destination / name
                shutil.copy2(work / name, target)
                target.chmod(0o755)


class NativeDistribution(Distribution):
    """Install the package at the wheel's platform-library root."""

    def has_ext_modules(self):
        return True


class BinaryWheel(bdist_wheel):
    """Tag the executable by platform, independent of the CPython ABI."""

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        _, _, platform = super().get_tag()
        return "py3", "none", platform


setup(distclass=NativeDistribution, cmdclass={"build_py": BuildFerre, "bdist_wheel": BinaryWheel})
