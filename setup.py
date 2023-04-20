#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages

author = 'James Robertson'

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


exec(open('cladeomatic/version.py').read())

setup(
    name='cladeomatic',
    include_package_data=True,
    version='0.0.1',
    python_requires='>=3.8.0,<4',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/cladeomaticc',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@phac-aspc.gc.ca',
    description=(
        'Clade-O-Matic: Automatic recognition of population structures based on canonical SNPs'),
    keywords='Genotyping, population structure, kmer',
    classifiers=classifiers,
    package_dir={'cladeomatic': 'cladeomatic'},
    package_data={
        "": ["*.txt", "*.fasta","*.html","*.gb"],
    },

    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'scipy',
        'six',
        'matplotlib',
        'argparse',
        'statistics',
        'plotly',
        'ete3',
        'jellyfish',
        'DendroPy',
        'pyahocorasick',
        'PyQt5',
        'ray',
        'deprecated'

    ],

    entry_points={
        'console_scripts': [
            'cladeomatic=cladeomatic.main:main',
        ],
    },
)