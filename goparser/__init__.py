import pkg_resources

from goparser.term import GOTerm
from goparser.annotation import GOAnnotation
from goparser.parser import GOParser

__version__ = pkg_resources.require('goparser')[0].version

__all__ = ['GOTerm','GOAnnotation','GOParser']
