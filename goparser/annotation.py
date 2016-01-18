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
from collections import Iterable

from goparser import GOTerm

class GOAnnotation(object):

    """Class representing an annotation of a gene with a GO term.

    For a list of annotation properties, see the
    `GAF 2.0 file format specification`__. 

    __ gafformat_

    Parameters
    ----------
    gene: str
        See :attr:`gene` attribute.
    term: `GOTerm` object
        See :attr:`term` attribute.
    evidence: str
        See :attr:`evidence` attribute.
    db_id: str, optional
        See :attr:`db_id` attribute.
    db_ref: list of str, optional
        See :attr:`db_ref` attribute.
    with_: list of str, optional
        See :attr:`with_` attribute.

    Attributes
    ----------
    gene: str
        The gene that is annotated (e.g., "MYOD1").
    term: `GOTerm` object
        The GO term that the gene is annotated with.
    evidence: str
        The three-letter evidence code of the annotation (e.g., "IDA").
    db_id: str, optional
        Database Object ID of the annotation.
    db_ref: list of str, optional
        DB:Reference of the annotation.
    with_: list of str, optional
        "With" information of the annotation.

    Methods
    -------
    get_gaf_format()
        Return the annotation as a tab-delimited string acccording to the
        `GAF 2.0 file format`__.

    get_pretty_format()
        Return a nicely formatted string representation of the GO annotation.


    __ gafformat_

    .. _gafformat: http://geneontology.org/page/go-annotation-file-gaf-format-20

    """

    _evidence_name = {\
            'EXP': 'experiment',\
            'IDA': 'direct assay',\
            'IPI': 'physical interaction',\
            'IMP': 'mutant phenotype',\
            'IGI': 'genetic interaction',\
            'IEP': 'expression pattern',\
            'ISS': 'sequence or structural similarity',\
            'ISO': 'sequence orthology',\
            'ISA': 'sequence alignment',\
            'ISM': 'sequence model',\
            'IGC': 'genomic context',\
            'IBA': 'biological aspect of ancestor',\
            'IBD': 'biological aspect of descendant',\
            'IKR': 'key residues',\
            'IRD': 'rapid divergence',\
            'RCA': 'reviewed computational analysis',\
            'TAS': 'traceable author statement',\
            'NAS': 'non-traceable author statement',\
            'IC' : 'inferred by curator',\
            'ND' : 'no biological data available',\
            'IEA': 'inferred from electronic annotation'\
            }
    """Mapping of the three-letter evidence codes to their full names.
    """

    _evidence_type = {\
            'EXP': 'experimental',\
            'IDA': 'experimental',\
            'IPI': 'experimental',\
            'IMP': 'experimental',\
            'IGI': 'experimental',\
            'IEP': 'experimental',\
            'ISS': 'computational',\
            'ISO': 'computational',\
            'ISA': 'computational',\
            'ISM': 'computational',\
            'IGC': 'computational',\
            'IBA': 'computational',\
            'IBD': 'computational',\
            'IKR': 'computational',\
            'IRD': 'computational',\
            'RCA': 'computational',\
            'TAS': 'literature',\
            'NAS': 'literature',\
            'IC' : 'curator',\
            'ND' : 'no_data',\
            'IEA': 'automatic'\
            }
    """Mapping of the three-letter evidence codes to their evidence types.
    """

    _evidence_type_short = {\
            'experimental': 'exp.',\
            'computational': 'comp.',\
            'literature': 'lit.',\
            'curator': 'cur.',\
            'no_data': 'n.d.',\
            'automatic': 'autom.'\
            }
    """Mapping of the evidence types to abbreviated forms.
    """

    #uniprot_pattern = re.compile("([A-Z][A-Z0-9]{5})(?:-(\d+))?")

    def __init__(self, gene, term, evidence, db_id = None,
                db_ref = None, with_ = None):

        assert isinstance(term, GOTerm)
        assert isinstance(gene, (str, unicode)) and gene != ''
        assert isinstance(evidence, (str, unicode)) and evidence != ''

        if db_id is not None:
            assert isinstance(db_id, (str, unicode)) and db_id != ''

        if db_ref is not None:
            assert isinstance(db_ref, Iterable)
        else:
            db_ref = []

        if with_ is not None:
            assert isinstance(with_, Iterable)

        self.gene = gene
        self.term = term
        self.evidence = evidence
        self.db_id = db_id
        self.db_ref = () if db_ref is None else tuple(db_ref)
        self.with_ = () if with_ is None else tuple(with_)

    def __repr__(self):
        return '<GOAnnotation (hash=%d)>' %(hash(self))

    def __str__(self):
        return '<GOAnnotation of gene "%s" with term "%s" (%s)>' \
                %(self.gene, self.term.name, self.term.id)

    def __eq__(self,other):
        if type(self) != type(other):
            return False
        elif self is other:
            return True
        else:
            return repr(self) == repr(other)

    def __hash__(self):
        data = []
        data.append(self.gene)
        data.append(self.term)
        data.append(self.evidence)
        data.append(self.db_id)
        data.append(self.db_ref)
        
        return hash(tuple(data))

    def get_gaf_format(self):
        """Return a GAF 2.0-compatible string representation of the annotation.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The formatted string.
        """
        sep = '\t'
        return sep.join([self.gene, self.db_ref, self.term.id, self.evidence,
                '|'.join(self.db_ref), '|'.join(self.with_)])

    #def get_pretty_format(self):
    #    """Returns a nicely formatted string with the annotation information.
    #
    #    Parameters
    #    ----------
    #    None
    #
    #    Returns
    #    -------
    #    str
    #        The formatted string.
    #    """
    #    pretty = "Annotation of gene '%s' with GO term '%s' (%s, reference: %s)'" \
    #            %(self.target,self.term.get_pretty_format(),self.evidence,'|'.join(self.db_ref))
    #    return pretty
