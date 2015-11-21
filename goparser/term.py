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

    """Class representing a GO term.

    This class is used by :func:`GOParser.parse_ontology` to store all parsed
    GO term data.

    Parameters
    ----------
    id_: str
        See ``id`` attribute.
    name: str
        See ``name`` attribute.
    domain: str
        See ``domain`` attribute.
    is_a: List of str
        See ``is_a`` attribute.
    part_of: List of str
        See ``part_of`` attribute.
    children:

    Attributes
    ----------
    id: str
        The ID of the GO term.
    name: str
        The name of the GO term.
    domain: str
        The domain of the GO term (e.g., "biological_process").
    is_a: set of str
        Set of GO term IDs that this GO term is a "subtype" of.
    part_of: set of str
        Set of GO term IDs that this GO term is a "part" of. 
    ancestors: set of str
        Set of GO term IDs that are "ancestors" of this GO term.
    children: set of str
        Set of GO term IDs that are "children" of this GO term.
    parts: set of str
        Set of GO term IDs that are "parts" of this GO term.
    descendants: set of str
        Set of GO terms IDs that are "descendants" of this GO term.

    Methods
    -------
    get_pretty_format(omit_acc=False, max_name_length=0, abbreviate=True)
        Returns a formatted version of the GO term name and ID.

    """

    _short_domain = {
        'biological_process': 'BP',
        'molecular_function': 'MF',
        'cellular_component': 'CC'
    }
    """Dictionary representing the abbreviations of the Gene Ontology domains.
    """

    _abbrev = [
        ('positive ','pos. '),
        ('negative ','neg. '),
        ('interferon-','IFN-'),
        ('proliferation','prolif.'),
        ('signaling','signal.')
    ]
    """List of tuples defining abbreviations to use in GO term names.
    """

    def __init__(self,id_,name,domain,is_a,part_of):

        self.id = id_ # unique identifier
        self.name = name
        self.domain = domain

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

    @staticmethod
    def id2acc(id_):
        """Converts a GO term ID to an accession number.

        Parameters
        ----------
        id_: str
            A GO term ID.

        Returns
        -------
        int
            The accession number corresponding to the GO term ID.

        """
        return int(id_[3:])

    @staticmethod
    def acc2id(self,acc):
        """Converts a GO term accession number to an ID.

        Parameters
        ----------
        acc: int
            A GO term accession number.

        Returns
        -------
        str
            The ID corresponding to the GO term accession number.
        """
        return 'GO:%07d' %(acc)

    @property
    def acc(self):
        return id2acc(self.id)

    @property
    def domain_short(self):
        return self._short_domain[self.domain]

    def get_pretty_format(self,omit_acc=False,max_name_length=0,abbreviate=True):
        """Returns a nicely formatted string with the GO term information.

        Parameters
        ----------
        omit_acc: bool, optional
            If set to True, don't include the GO term ID.
        max_name_length: int, optional
            If set, the formatted string (excluding the ID) will be truncated
            so that its total length does not exceed this value.
        abbreviate: bool, optional
            If set to False, do not use abberviations (see ``_abbrev``) to
            shorten the GO term name.

        Returns
        -------
        str
            The formatted string.
        """
        name = self.name
        if abbreviate:
            for abb in self._abbrev:
                name = re.sub(abb[0],abb[1],name)
        if max_name_length >= 3 and len(name) > max_name_length:
            name = name[:(max_name_length-3)] + '...'
        if omit_acc: return "%s: %s" %(self.domain_short, name)
        else: return "%s: %s (%s)" %(self.domain_short, name, self.id)

    def get_tuple(self):
        """Returns a 4-tuple containing the GO term information.

        Parameters
        ----------
        None

        Returns
        -------
        tuple (of length 4)
            A tuple with elements consisting of the GO term ID, the string
            "GO", the shortened domain (see ``_short_domain``), and the GO term
            name.
        """
        return (self.id,'GO',self.domain_short,self.name)


