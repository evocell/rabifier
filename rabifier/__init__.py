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

import os
import json
import logging
import tempfile


__title__ = 'rabifier'
__version__ = '2.0.2'
__author__ = 'Jaroslaw Surkont'
__license__ = 'GPLv3'
__copyright__ = 'Copyright 2015 Jaroslaw Surkont'

# load configuration
with open(os.path.join(os.path.dirname(__file__), 'config.json')) as fin:
    config = json.load(fin)
config['tmp'] = tempfile.gettempdir()
config['seed']['path'] = os.path.join(os.path.dirname(__file__), 'data')

# Configure logging
#logging.basicConfig(level=logging.DEBUG)
logging.getLogger(__name__)
