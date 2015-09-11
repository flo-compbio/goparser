# Copyright (c) 2015 Florian Wagner
#
# This file is part of GOParser.
#
# GOParser is free software: you can redistribute it and/or modify
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
	#abbrev = [('regulation','reg.'),('positive','pos.'),('negative','neg.'),('interferon-','IFN-')]

	def __init__(self,id_,name,namespace,is_a,part_of,definition=None):
		self.id = id_
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
		#print self.namespace
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


class GOAnnotation(object):

	"""
	Class representing a GO annotation (i.e., an association of a GO term with a gene).
	"""

	#uniprot_pattern = re.compile("([A-Z][A-Z0-9]{5})(?:-(\d+))?")
	#short_ns = {'biological_process': 'BP', 'molecular_function': 'MF', 'cellular_component': 'CC'}

	def __init__(self,target,term,evidence,db_id=None,db_ref=[],with_=[]):
		assert target is not None and target != ''
		assert evidence is not None and evidence != ''
		assert type(term) == GOTerm
		self.target = target # a gene name
		self.evidence = evidence # GO evidence code
		self.term = term # a GOTerm object
		#self.pubmed_id = pubmed_id # PubMed ID, optional
		#self.uniprot = uniprot # Uniprot identifier, optional
		self.db_id = db_id
		self.db_ref = db_ref[:]
		self.with_ = with_[:]

	evidence_name = {\
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

	evidence_type = {\
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

	evidence_type_short = {\
			'experimental': 'exp.',\
			'computational': 'comp.',\
			'literature': 'lit.',\
			'curator': 'cur.',\
			'no_data': 'n.d.',\
			'automatic': 'autom.'\
			}

	def __repr__(self):
		return "<GOAnnotation %s>" %(', '.join([self.target,self.db_id,self.evidence,'|'.join(self.db_ref),'|'.join(self.with_),repr(self.term)]))

	def __str__(self):
		#return "<GOAnnotation: %s>" %(self.get_formatted(sep=','))
		return "<GOAnnotation: %s>" %(self.get_pretty_format())

	def __eq__(self,other):
		if type(self) != type(other):
			return False

		elif self is other:
			return True

		elif self.target == other.target \
				and self.db_id == other.db_id \
				and self.term == other.term \
				and self.evidence == other.evidence \
				and self.db_ref == other.db_ref \
				and self.with_ == other.with_:
			return True

		else:
			return False

	def __hash__(self):
		return hash(repr(self))

	def get_formatted(self,sep='\t'):
		return sep.join([self.target,self.db_ref,self.term.id,self.evidence,'|'.join(self.db_ref),'|'.join(self.with_)])

	def get_pretty_format(self):
		pretty = "Annotation of gene '%s' with GO term '%s' (%s, reference: %s)'" \
				%(self.target,self.term.get_pretty_format(),self.evidence,'|'.join(self.db_ref))
		return pretty

