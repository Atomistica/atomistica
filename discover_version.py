# Copyright © 2019 Lars Pastewka
#
# µSpectre is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Lesser Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# µSpectre is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with µSpectre; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#
# Additional permission under GNU GPL version 3 section 7
#
# If you modify this Program, or any covered work, by linking or combining it
# with proprietary FFT implementations or numerical libraries, containing parts
# covered by the terms of those libraries' licenses, the licensors of this
# Program grant you additional permission to convey the resulting work.

#
# This is the most minimal-idiotic way of discovering the version that I
# could come up with. It deals with the following issues:
# * If we are installed, we can get the version from package metadata,
#   either via importlib.metadata or from pkg_resources. This also holds for
#   wheels that contain the metadata. We are good! Yay!
# * If we are not installed, there are two options:
#   - We are working within the source git repository. Then
#        git describe --tags --always
#     yields a reasonable version descriptor, but that is unfortunately not
#     PEP 440 compliant (see https://peps.python.org/pep-0440/). We need to
#     mangle the version string to yield something compatible.
# - If we are not in a source directory, check if muGrid's version.cc has
#   been created and extract version from that file.
# - If we install from a source tarball, all version information may be lost.
#   Fortunately, Meson uses git archive to create the source tarball, which
#   replaces certain tags with commit information. Unfortunately, what this
#   yields is different from `git describe` - in particular, it only yields the
#   tag (which contains the version information) if we are *exactly* on the
#   tag commit. (`git describe` tells us the distance from the latest tag.) We
#   need to extract the version information from the string provided, but if
#   we are not on the tag we can only return a bogus version (here 0.0.0.0).
#   It works for releases, but not for a tarball generated from a random
#   commit. I am not happy and open for suggestions.
#

import re
import subprocess
import sys

class CannotDiscoverVersion(Exception):
    pass


def get_version_from_pkg_info():
    """
    Discover version from PKG-INFO file. This is never dirty and we cannot
    obtain the git hash. This is typically no problem since the version
    uniquely identifies the commit.
    """
    try:
        with open('PKG-INFO', 'r') as f:
            l = f.readline()
            while l:
                if l.startswith('Version:'):
                    return False, l[8:].strip(), '(PKG-INFO)'
                l = f.readline()
    except FileNotFoundError:
        raise CannotDiscoverVersion("File 'PKG-INFO' does not exist.")
    raise CannotDiscoverVersion("No line starting with 'Version:' in 'PKG-INFO'.")


def get_version_from_git():
    """
    Discover muSpectre version from git repository. We get all information,
    including whether the version is dirty and the git hash.
    """
    git_describe = subprocess.run(
        ['git', 'describe', '--tags', '--dirty', '--always'],
        stdout=subprocess.PIPE)
    if git_describe.returncode != 0:
        raise CannotDiscoverVersion('git execution failed')
    version = git_describe.stdout.decode('latin-1').strip()
    git_show = subprocess.run(
        ['git', 'show', '-s', '--format=%H'],
        stdout=subprocess.PIPE)
    if git_show.returncode != 0:
        raise CannotDiscoverVersion('git execution failed')
    git_hash = git_show.stdout.decode('latin-1').strip()

    dirty = version.endswith('-dirty')

    # Make version PEP 440 compliant
    if dirty:
        version = version.replace('-dirty', '')
    version = version.replace('-', '+', 1)
    if dirty:
        version += '-dirty'

    return dirty, version, git_hash


def get_version_from_cc(fn):
    """
    Extract version from a version.cc file (that is typically created by the)
    build system. This will also yield dirty state and git hash, if present
    in the cc file.
    """
    text = open(fn, 'r').read()
    dirty = bool(re.search('constexpr bool git_dirty{(true|false)};',
                           text).group(1))
    version = re.search(
        r'constexpr char git_describe[]{"([A-Za-z0-9_.+-]*)"};', text
    ).group(1)
    git_hash = re.search(r'constexpr char git_hash[]{"([A-Za-z0-9_.-]*)"};',
                         text).group(1)

    return dirty, version, git_hash


#
# Discover version and write version.cc if version comes from git
#

try:
    try:
        # First try is to look inside PKG-INFO file.
        dirty, version, git_hash = get_version_from_pkg_info()
    except CannotDiscoverVersion:
        dirty, version, git_hash = get_version_from_git()
except CannotDiscoverVersion:
    # Detection via PKG-INFO and git failed. Get version from version.cc file.
    dirty, version, git_hash = \
        get_version_from_cc('src/libmugrid/version.cc')

#
# Print version to screen
#

if len(sys.argv) > 1 and '--full' in sys.argv:
    # True and false need to be C++ compatible
    print('true' if dirty else 'false', version, git_hash)
else:
    print(version)
