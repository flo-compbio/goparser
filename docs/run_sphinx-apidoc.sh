#!/bin/bash

sphinx-apidoc -e -o source/api ../goparser
rm source/api/modules.rst
rm source/api/goparser.rst
