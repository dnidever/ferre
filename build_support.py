"""Platform settings shared by FERRE's source compiler and wheel tags."""
import os
import re
import subprocess
import sys
from pathlib import Path


def macos_options(flags, link_flags):
    """Return explicit SDK/target flags and the corresponding wheel target.

    A Python built for an old macOS release can report a deployment tag older
    than the system/compiler being used. Use one target for compilation and the
    wheel name; dependency auditing by delocate still determines whether all
    bundled runtime libraries actually support that target.
    """
    if sys.platform != 'darwin':
        return list(flags), list(link_flags), None
    flags, link_flags = list(flags), list(link_flags)
    explicit = [f.split('=', 1)[1] for f in flags + link_flags
                if f.startswith('-mmacosx-version-min=')]
    target = explicit[-1] if explicit else os.environ.get('MACOSX_DEPLOYMENT_TARGET')
    if not target:
        version = subprocess.check_output(['sw_vers', '-productVersion'], text=True).strip()
        # Big Sur is understood by older GNU Fortran compilers as well as
        # modern ones. CI explicitly selects the audited release target.
        major, *minor = version.split('.')
        target = '11.0' if int(major) >= 11 else major + '.' + minor[0]
    if not re.fullmatch(r'\d+\.\d+(?:\.\d+)?', target):
        raise RuntimeError('MACOSX_DEPLOYMENT_TARGET must be a macOS version such as 13.0')
    if not explicit:
        flags.append('-mmacosx-version-min=' + target)
    if not any('syslibroot' in f or f.startswith(('-isysroot', '--sysroot'))
               for f in flags + link_flags):
        sdk = os.environ.get('SDKROOT')
        if not sdk:
            try:
                sdk = subprocess.check_output(['xcrun', '--show-sdk-path'], text=True).strip()
            except (OSError, subprocess.CalledProcessError) as error:
                raise RuntimeError('Cannot locate the macOS SDK. Install Xcode Command Line Tools '
                                   'or set SDKROOT/LDFLAGS explicitly.') from error
        if not Path(sdk).is_dir():
            raise RuntimeError('The selected macOS SDK directory does not exist: ' + sdk)
        link_flags.append('-Wl,-syslibroot,' + sdk)
    return flags, link_flags, '.'.join(target.split('.')[:2])
