# Copyright (c) 2015, 2016 Florian Wagner
#
# This file is part of GOparser.
#
# GOparser is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import pkg_resources

from goparser.term import GOTerm
from goparser.annotation import GOAnnotation
from goparser.parser import GOParser

__version__ = pkg_resources.require('goparser')[0].version

__all__ = ['GOTerm', 'GOAnnotation', 'GOParser']
