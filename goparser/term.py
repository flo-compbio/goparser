# Copyright (c) 2015 Florian Wagner
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

import re

class GOTerm(object):

    """
    Class representing a GO term.
    """

    short_ns = {'biological_process': 'BP', 'molecular_function': 'MF', 'cellular_component': 'CC'}
    abbrev = [('positive ','pos. '),('negative ','neg. '),('interferon-','IFN-'),('proliferation','prolif.'),('signaling','signal.')]

    def __init__(self,id_,name,namespace,is_a,part_of,definition=None):

        self.id = id_ # unique identifier
        self.name = name
        self.namespace = namespace
        self.definition = definition

        # to store immediate parents/wholes
        self.is_a = is_a.copy()
        self.part_of = part_of.copy()

        # to store immediate children/parts
        self.children = set()
        self.parts = set()

        # to store all descendants/ancestors
        self.descendants = None
        self.ancestors = None


    def __repr__(self):
        return "<GOTerm %s>" %(self.id)
    def __str__(self):
        return "<GOTerm: %s>" %(self.get_pretty_format())

    def __eq__(self,other):
        if type(self) != type(other):
            return False

        elif self is other:
            return True
            
        elif self.id == other.id:
            return True

        else:
            return False

    def __hash__(self):
        return hash(repr(self))

    def get_namespace_short(self):
        return self.short_ns[self.namespace]

    def get_id(self):
        return self.id

    def get_acc(self): # gets the accession number as integer
        return int(self.id[3:])

    def get_pretty_format(self,omit_acc=False,max_name_length=0,abbreviate=True):
        name = self.name
        if abbreviate:
            for abb in self.abbrev:
                name = re.sub(abb[0],abb[1],name)
        if max_name_length >= 3 and len(name) > max_name_length:
            name = name[:(max_name_length-3)] + '...'
        if omit_acc: return "%s: %s" %(self.short_ns[self.namespace], name)
        #else: return "%s: %s (GO:%07d)" %(self.short_ns[self.namespace], self.name, self.acc)
        else: return "%s: %s (%s)" %(self.short_ns[self.namespace], name, self.id)

    def get_tuple(self):
        return (self.id,'GO',self.get_namespace_short(),self.name)


