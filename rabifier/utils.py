#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of Rabifier.
#
# Rabifier is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Rabifier is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rabifier.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import os
import subprocess
import shutil
import logging


logger = logging.getLogger(__name__)


class Pathfinder(object):

    def __init__(self, include_sys, *paths):
        """

        :param include_sys: include system paths
        :type include_sys: bool
        :param paths: paths to additional directories for searching
        """

        self.paths = set()
        for path in paths:
            self.add_path(path)
        if include_sys:
            for path in os.environ['PATH'].split(os.pathsep):
                self.add_path(path)

    def add_path(self, path):
        """ Adds a new path to the list of searchable paths

        :param path: new path
        """

        if os.path.exists(path):
            self.paths.add(path)
            return path
        else:
            #logger.debug('Path {} doesn\'t exist'.format(path))
            return None

    def __getitem__(self, item):
        return self.get(item)

    def get(self, name):
        """ Looks for a name in the path.

        :param name: file name
        :return: path to the file
        """

        for d in self.paths:
            if os.path.exists(d) and name in os.listdir(d):
                return os.path.join(d, name)
        logger.debug('File not found {}'.format(name))
        return None

    def exists(self, name):
        """ Checks if a file matching the name exist in the path

        :param name: file name
        :return: True if file exists else False
        """

        return False if self.get(name) is None else True


def run_cmd(cmd, out=os.path.devnull, err=os.path.devnull):
    """Runs an external command

    :param list cmd: Command to run.
    :param str out: Output file
    :param str err: Error file
    :raises: RuntimeError
    """

    logger.debug(' '.join(cmd))
    with open(out, 'w') as hout:
        proc = subprocess.Popen(cmd, stdout=hout, stderr=subprocess.PIPE)
        err_msg = proc.communicate()[1].decode()
    with open(err, 'w') as herr:
        herr.write(str(err_msg))
    msg = '({}) {}'.format(' '.join(cmd), err_msg)
    if proc.returncode != 0:
        logger.error(msg)
        raise RuntimeError(msg)


def run_cmd_if_file_missing(cmd, fname, out=os.path.devnull, err=os.path.devnull):
    """Runs an external command if file is absent.

    :param list cmd: Command to run.
    :param str fname: Path to the file, which existence is being checked.
    :param str out: Output file
    :param str err: Error file
    :return: True if cmd was executed, False otherwise
    :rtype: boolean
    """

    if fname is None or not os.path.exists(fname):
        run_cmd(cmd, out, err)
        return True
    else:
        return False


def merge_files(sources, destination):
    """Copy content of multiple files into a single file.

    :param list(str) sources: source file names (paths)
    :param str destination: destination file name (path)
    :return:
    """

    with open(destination, 'w') as hout:
        for f in sources:
            if os.path.exists(f):
                with open(f) as hin:
                    shutil.copyfileobj(hin, hout)
            else:
                logger.warning('File is missing: {}'.format(f))


