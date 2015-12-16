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

"""Module containing the GOParser class.

"""

import csv
import gzip
import re
import sys
import logging
import bisect
import cPickle as pickle
from collections import Counter

from genometools import misc
from goparser import GOTerm, GOAnnotation

logger = logging.getLogger(__name__)

class GOParser(object):
    """ A class for accessing Gene Ontology (GO) term and annotation data.

    This class provides functions for parsing text files describing the Gene
    Ontology and GO annotations, for accessing information about specific GO
    terms, as well as for querying the data for associations between genes and
    GO terms.

    Parameters
    ----------
    quiet: bool, optional
        If True, only warnings and errors will be reported.
    verbose: bool, optional
        If True, enable verbose logging (i.e., including debug messages).
        If ``quiet`` is set to True, the value of this parameter is ignored.

    Attributes
    ----------
    terms: dict [str:GOTerm]
        A mapping of GO term IDs to `GOTerm` objects, each representing a
        single GO term. Populated by the member function `parse_ontology`.
    genes: set of str
        A set of all "valid" gene names. Populated by the member function
        `parse_annotations`. Typically, this is the set of all protein-coding
        genes of a particular species. GOparser ignores all annotations
        for genes that are not in this set.
    annotations: list of GOAnnotation objects
        A list of `GOAnnotation` objects, each representing a single GO
        annotation. Populated by the member function `parse_annotations`.
    term_annotations: dict [str:list of GOAnnotation objects]
        A mapping of GO term IDs to lists of `GOAnnotation` objects, with each
        list representing all annotations that use a particular GO term.
    gene_annotations: dict [str:list of GOAnnotation objects]
        A mapping of gene symbols to lists of `GOAnnotation` objects, with each
        list representing all annotations of a particular gene.

    Methods
    -------
    get_term_by_id(id_)
        Return the term with the given term ID as a `GOTerm` object.
    get_term_by_name(name)
        Return the term with the given name as a `GOTerm` object.
    get_gene_goterms(gene, ancestors=False)
        Return all GO terms that the given gene is annotated with.
        If ``ancestors`` is set to True, also return all ancestor GO terms
        of those terms.
    get_goterm_genes(id_, descendants=True)
        Return all genes annotated with the GO term corresponding to the given
        GO term ID. If ``descendants`` is set to True, also
        return genes annotated with any descendant GO term of this term. Since
        annotations should be propagated down to descendant terms, this is the
        default behavior.
    save(ofn, compress=False)
        Stores the GOParser object as a `pickle` file. If ``compress`` is set
        to True, the object is stored as a gzip'ed pickle file.
    load(fn)
        Loads the GOParser object from a `pickle` file. Gzip compression is
        detected automatically.
        
    Notes
    -----
    The typical workflow for reading the GO annotations for a specific species
    looks as follows:

    - Step 1) Extract a list of protein-coding genes using the script
      ``extract_protein_coding_genes.py`` from the `genometools` package.
      (See the `GenomeTools documentation`__.)

    __ extract_genes_

    - Step 2) Use the `parse_ontology` member function to parse the
      ``go-basic.obo`` file, containing the Gene Ontology.
      (This file can be `downloaded`__ from the website of the Gene Ontology
      Consortium.)

    __ download_go_

    - Step 3) Use the `parse_annotations` member function to parse a gene
      association file (GAF), containing annotations of genes with GO terms for
      a specific species. The list of protein-coding genes generated in Step 1)
      is used to only parse annotations for protein-coding genes.
      A species-specific file can be `downloaded`__ from the ftp server of the
      UniProt-GOA database.

    __ download_gaf_

    Afterwards, the member functions and `get_term_by_id` and
    `get_term_by_name` can be used to obtain GOTerm objects containing
    information about individual GO terms. The member function
    `get_gene_goterms` can be used to obtain a list of all GO terms a particular
    gene is annoatated with, and the member function `get_goterm_genes` can be
    used to obtain a list of all genes annotated with a particular GO term.

    .. _extract_genes: https://genometools.readthedocs.org/en/latest/scripts.html#extract-protein-coding-genes-py
    .. _download_go: http://geneontology.org/page/download-ontology
    .. _download_gaf: http://www.ebi.ac.uk/GOA/downloads

    Examples
    --------
    The following example assumes that the Gene Ontology OBO file and the
    UniProt-GOA gene association files have been downloaded, and that a list of
    protein-coding genes named "protein_coding_genes_human.tsv" has been
    generated using the genometools Python package.

    >>> from goparser import GOParser
    >>> G = GOParser()
    >>> GOParser.parse_ontology('go-basic.obo')
    >>> GOParser.parse_annotations('gene_association.goa_human.gz','protein_coding_genes_human.tsv')
    >>> print GOParser.get_gene_goterms('MYC')

    """

    def __init__(self, quiet = False, verbose = False):
        self.terms = {}
        self.annotations = []
        self.term_annotations = {}
        self.gene_annotations = {}

        self._syn2id = {}
        self._alt_id = {}
        self._name2id = {}
        self._flattened = False

    def save(self, ofn, compress = False):
        """Serialize the current GOParser object and store it in a pickle file.

        Parameters
        ----------
        ofn: str
            Path of the output file.
        compress: bool, optional
            Whether to compress the file using gzip.

        Returns
        -------
        None

        Notes
        -----
        Compression with gzip is significantly slower than storing the file
        in uncompressed form.
        """
        logger.info('Saving pickle...')
        if compress:
            with gzip.open(ofn, 'wb') as ofh:
                pickle.dump(self, ofh, pickle.HIGHEST_PROTOCOL)
        else:
            with open(ofn, 'wb') as ofh:
                pickle.dump(self, ofh, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(fn):
        """Load a GOParser object from a pickle file.

        The function automatically detects whether the file is compressed
        with gzip.

        Parameters
        ----------
        fn: str
            Path of the pickle file.

        Returns
        -------
        GOParser
            The GOParser object stored in the pickle file.
        """
        P = None
        with misc.open_plain_or_gzip(fn, 'rb') as fh:
            P = pickle.load(fh)
        return P

    def get_term_by_id(self, id_):
        """Get the GO term corresponding to the given GO term ID.

        Parameters
        ----------
        id_: str
            A GO term ID.

        Returns
        -------
        GOTerm
            The GO term corresponding to the given ID.
        """
        return self.terms[id_]

    def get_term_by_acc(self, acc):
        """Get the GO term corresponding to the given GO term accession number.

        Parameters
        ----------
        acc: int
            The GO term accession number.

        Returns
        -------
        GOTerm
            The GO term corresponding to the given accession number.
        """
        return self.terms[GOParser.acc2id(acc)]

    def get_term_by_name(self, name):
        """Get the GO term with the given GO term name.

        If the given name is not associated with any GO term, the function will
        search for it among synonyms.

        Parameters
        ----------
        name: str
            The name of the GO term.

        Returns
        -------
        GOTerm
            The GO term with the given name.

        Raises
        ------
        ValueError
            If the given name is found neither among the GO term names, nor
            among synonyms.
        """
        term = None
        func_name = 'get_term_by_name'
        try:
            term = self.terms[self._name2id[name]]
        except KeyError:
            try:
                term = self.terms[self._syn2id[name]]
            except KeyError:
                pass
            else:
                logger.warning('%s: GO term name "%s" is a synonym for "%s".',
                        func_name, name, term.name)

        if term is None:
            raise ValueError('%s : GO term name "%s" not found!'
                    %(func_name, name))

        return term

    def clear_data(self):
        """Clear both ontology and annotation data.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.clear_annotation_data()
        self.terms = {}
        self._alt_id = {}
        self._syn2id = {}
        self._name2id = {}
        self._flattened = False
 
    def clear_annotation_data(self):
        """Clear annotation data.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.genes = set()
        self.annotations = []
        self.term_annotations = {}
        self.gene_annotations = {}

    def parse_ontology(self, fn, flatten = True, part_of_cc_only = False):
        """ Parse an OBO file and store GO term information.

        This function needs to be called before `parse_annotations`, in order
        to read in the Gene Ontology terms and structure.

        Parameters
        ----------
        fn: str
            Path of the OBO file.
        flatten: bool, optional
            If set to False, do not generate a list of all ancestors and
            descendants for each GO term. Warning: Without flattining,
            GOparser cannot propagate GO annotations properly.
        part_of_cc_only: bool, optional
            Legacy parameter for backwards compatibility. If set to True,
            ignore ``part_of`` relations outside the ``celluclar_component``
            domain.

        Notes
        -----
        The function erases all previously parsed data.
        The function requires the OBO file to end with a line break.
        """
        self.clear_data() # clear all old data

        with open(fn) as fh:
            n = 0
            while True:
                try:
                    nextline = fh.next()
                except StopIteration:
                    break
                if nextline == '[Term]\n':
                    n+=1
                    id_ = fh.next()[4:-1]
                    #acc = get_acc(id_)
                    name = fh.next()[6:-1]
                    self._name2id[name] = id_
                    domain = fh.next()[11:-1]
                    is_a = set()
                    part_of = set()
                    l = fh.next()
                    while l != '\n':
                        if l.startswith('alt_id:'):
                            self._alt_id[l[8:-1]] = id_
                        elif l.startswith('is_a:'):
                            is_a.add(l[6:16])
                        elif l.startswith('synonym:'):
                            idx = l[10:].index('"')
                            if l[(10+idx+2):].startswith("EXACT"):
                                s = l[10:(10+idx)]
                                self._syn2id[s] = id_
                        elif l.startswith('relationship: part_of'):
                            if part_of_cc_only:
                                if domain == 'cellular_component':
                                    part_of.add(l[22:32])
                            else:
                                part_of.add(l[22:32])
                        l = fh.next()
                    self.terms[id_] = GOTerm(id_, name, domain, is_a, part_of)

        logger.info('Parsed %d GO term definitions.', n)

        # store children and parts
        logger.info('Adding child and part relationships...')
        for id_,term in self.terms.iteritems():
            for parent in term.is_a:
                self.terms[parent].children.add(id_)
            for whole in term.part_of:
                self.terms[whole].parts.add(id_)

        if flatten:
            logger.info('Flattening ancestors...')
            self._flatten_ancestors()
            logger.info('Flattening descendants...')
            self._flatten_descendants()
            self._flattened = True

    def _flatten_ancestors(self, include_part_of = True):
        """Determines and stores all ancestors of each GO term.

        Parameters
        ----------
        include_part_of: bool, optional
            Whether to include ``part_of`` relations in determining
            ancestors.

        Returns
        -------
        None
        """

        def get_all_ancestors(term):
            ancestors = set()
            for id_ in term.is_a:
                ancestors.add(id_)
                ancestors.update(get_all_ancestors(self.terms[id_]))
            if include_part_of:
                for id_ in term.part_of:
                    ancestors.add(id_)
                    ancestors.update(get_all_ancestors(self.terms[id_]))
            return ancestors

        for term in self.terms.itervalues():
            term.ancestors = get_all_ancestors(term)

    def _flatten_descendants(self, include_parts = True):
        """Determines and stores all descendants of each GO term.

        Parameters
        ----------
        include_parts: bool, optional
            Whether to include ``part_of`` relations in determining
            descendants.

        Returns
        -------
        None
        """

        def get_all_descendants(term):
            descendants = set()
            for id_ in term.children:
                descendants.add(id_)
                descendants.update(get_all_descendants(self.terms[id_]))
            if include_parts:
                for id_ in term.parts:
                    descendants.add(id_)
                    descendants.update(get_all_descendants(self.terms[id_]))
            return descendants

        for term in self.terms.itervalues():
            term.descendants = get_all_descendants(term)

    def parse_annotations(self, annotation_file, gene_file,
            db_sel = 'UniProtKB',
            select_evidence = None, exclude_evidence = None,
            exclude_ref = None, strip_species = False, ignore_case = False):
        """Parse a GAF file and store annotations as GOAnnotation objects.

        Parameters
        ----------
        annotation_file: str
            Path of the annotation file.
        gene_file: str
            Path of the gene file. (This file can be generated using
            extract_protein_coding_genes.py from `genometools`.)
        db_sel: str, optional
            Select only annotations with this ``DB`` (column 1) value.
            If empty, disable filtering based on the ``DB`` value.
        select_evidence: list of str, optional
            Only include annotations with the given evidence codes.
            It not specified, allow all evidence codes, except for those listed
            in ``exclude_evidence``.
        exclude_evidence: list of str, optional
            Exclude all annotations with any of the given evidence codes.
            If ``select_evidence`` is specified, this parameter is ignored.
            If not specified, allow all evidence codes.
        exclude_ref: list of str, optional
            Exclude all annotations with the given DB:reference (column 6).
            Example: ``["PMID:2676709"]``. Note: This filter is currently
            ignored if an annotation has more than one reference.
        strip_species: bool, optional
            Undocumented.
        ignore_case: bool, optional
            Undocumented.

        Returns
        -------
        None
        """

        if not self.terms:
            raise ValueError('You need to first parse an OBO file!')

        if select_evidence is None:
            select_evidence = []

        if exclude_evidence is None:
            exclude_evidence = []

        if exclude_ref is None:
            exclude_ref = []

        # always overwrite all previously parsed annotations
        self.clear_annotation_data()

        # read genes
        genes = set()
        genes_upper = {}
        with open(gene_file) as fh:
            reader = csv.reader(fh,dialect='excel-tab')
            for l in reader:
                genes.add(l[0])
                if ignore_case:
                    genes_upper[l[0].upper()] = l[0]
        self.genes = genes # store the list of genes for later use
        logger.info('Read %d genes.', len(genes))

        # read annotations
        self.term_annotations = dict((id_,[]) for id_ in self.terms)
        self.gene_annotations = dict((g,[]) for g in self.genes)
        gene_terms = dict((g,set()) for g in self.genes) # used for statistics

        isoform_pattern = re.compile(r"UniProtKB:([A-Z][0-9A-Z]{5}-\d+)")
        gene_pattern = re.compile(r"[a-zA-Z0-9]+\.\d+$")
        pmid_pattern = re.compile(r"(?:PMID:\d+|DOI:[^\s]+)")
        uniprot_pattern = re.compile(r"UniProtKB:([A-Z][0-9A-Z]{5}(?:-\d+)?)")

        unknown_gene_names = Counter()
        unknown_gene_annotations = 0

        unknown_term_ids = Counter()
        unknown_term_annotations = 0

        # Parsing!
        logger.info('Parsing annotations...')
        n = 0
        excluded_evidence_annotations = 0
        excluded_reference_annotations = 0
        valid_annotations = 0
        with misc.smart_open(annotation_file, try_gzip = True) as fh:
            reader = csv.reader(fh,dialect='excel-tab')
            for i,l in enumerate(reader):
                target = None

                if not l: continue
                if ((not db_sel) or l[0] == db_sel) and l[3] != 'NOT':
                    n+=1

                    # test if evidence code is excluded
                    if (select_evidence and l[6] not in select_evidence) \
                            or l[6] in exclude_evidence:
                        excluded_evidence_annotations += 1
                        continue

                    # test if reference is excluded
                    db_ref = []
                    if l[5]:
                        db_ref = l[5].split('|')
                        if len(db_ref) == 1 and db_ref[0] in exclude_ref:
                            excluded_reference_annotations += 1
                            continue
                            

                    # determine target
                    if not l[2]:
                        raise Exception('Missing target gene in line %d:\n%s'
                                %(i+1, '\t'.join(l)))

                    gene = l[2]
                    db = l[0]
                    db_id = l[1]
                    if strip_species:
                        try:
                            gene = gene[:gene.rindex('_')]
                        except ValueError as e:
                            pass

                    term_id = l[4]
                    evidence = l[6]

                    invalid = False

                    if (ignore_case and gene.upper() not in genes_upper) \
                            or ((not ignore_case) and gene not in self.genes):
                        unknown_gene_annotations += 1
                        unknown_gene_names[l[2]] += 1
                        invalid = True

                    if term_id not in self.terms:
                        unknown_term_annotations += 1
                        unknown_term_ids[term_id] += 1
                        invalid = True

                    if not invalid:
                
                        valid_annotations += 1

                        # if ignore_case, convert gene to "original" name
                        if ignore_case:
                            gene = genes_upper[gene.upper()]

                        term = self.terms[term_id]

                        # parse secondary information
                        # (associated UniProt and PubMed entries)
                        #pmid = pmid_pattern.search(l[5])
                        #if pmid is not None: pmid = pmid.group(0)
                        #uniprot = uniprot_pattern.search(l[7])
                        #if uniprot is not None: uniprot = uniprot.group(1)
                        with_ = []
                        if l[7]:
                            with_ = l[7].split('|')

                        # generate annotation
                        ann = GOAnnotation(target = gene, term = term,
                                evidence = evidence, db_id = db_id,
                                db_ref = db_ref, with_ = with_)

                        # add annotation to global list
                        self.annotations.append(ann)

                        # add annotation under term ID
                        self.term_annotations[term_id].append(ann)

                        # add annotation under gene
                        self.gene_annotations[gene].append(ann)
                        gene_terms[gene].add(term_id)

        # output some statistics
        if n > 0:
            logger.info('Parsed %d positive GO annotations ' +
                    '(%d = %.1f%% excluded based on evidence type).', \
                    n, excluded_evidence_annotations,
                    100*(excluded_evidence_annotations / float(n)))

        if unknown_gene_annotations > 0:
            logger.warning('Warning: %d annotations with %d unkonwn gene ' +
                    'names.',
                    unknown_gene_annotations, len(unknown_gene_names))

        if unknown_term_annotations > 0:
            logger.warning('Warning: %d annotations with %d unkonwn term IDs.',
                    unknown_term_annotations, len(unknown_term_ids))

        logger.info('Found a total of %d valid annotations.',
                valid_annotations)

        logger.info('%d unique Gene-Term associations.',
                sum(len(gene_terms[g]) for g in genes))

    def get_gene_goterms(self, gene, ancestors = False):
        """Return all GO terms a particular gene is annotated with.

        Parameters
        ----------
        gene: str
            The gene symbol of the gene.
        ancestors: bool, optional
            If set to True, also return all ancestor GO terms.


        Returns
        -------
        set of GOTerm objects
            The set of GO terms the gene is annotated with.

        Notes
        -----
        If a gene is annotated with a particular GO term, it can also be
        considered annotated with all ancestors of that GO term.
        """
        annotations = self.gene_annotations[gene]
        terms = set(ann.term for ann in annotations)

        if ancestors:
            assert self._flattened
            ancestor_terms = set()
            for t in terms:
                ancestor_terms.update(self.terms[id_] for id_ in t.ancestors)
            terms |= ancestor_terms

        return terms

    def get_goterm_genes(self, id_, descendants = True):
        """Return all genes that are annotated with a particular GO term.

        Parameters
        ----------
        id_: str
            GO term ID of the GO term.
        descendants: bool, optional
            If set to False, only return genes that are directly annotated with
            the specified GO term. By default, also genes annotated with any
            descendant term are returned.

        Returns
        -------
        Notes
        """

        # determine which terms to include
        main_term = self.terms[id_]
        check_terms = set([main_term])

        if descendants:
            assert self._flattened
            check_terms.update([self.terms[id_]
                    for id_ in main_term.descendants])

        # get annotations of all included terms
        genes = set()
        for term in check_terms:
            genes.update(ann.target for ann in self.term_annotations[term.id])

        return genes
