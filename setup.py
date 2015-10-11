from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
description = 'GOparser - A Python framework for working with gene ontology (GO) terms and annotations'
version = '1.1rc2'

long_description = ''
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='goparser',

    version=version,

    description=description,
    long_description=long_description,

    # homepage.
    url='https://github.com/flo-compbio/goparser',

    author='Florian Wagner',
    author_email='florian.wagner@duke.edu',

    license='GPLv3',

    # see https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='gene ontology biology bioinformatics',

    #packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    packages=['goparser'],

    install_requires=['genometools','sphinx','sphinx_rtd_theme'],

	# development dependencies
    #extras_require={
    #},

	# data
    #package_data={
    #},

	# data outside package
    #data_files=[('my_data', ['data/data_file'])],

	# executable scripts
    #entry_points={
    #    'console_scripts': []
    #},
)
