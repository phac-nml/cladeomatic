#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages
from cladeomatic.version import __version__
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
    version=__version__,
    python_requires='>=3.9.0,<4',
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
        'numpy>=1.23.0',
        'pandas>=2.0.1',
        'biopython>=1.81',
        'six>=1.16.0',
        'matplotlib>=3.7.1',
        'statistics>=1.0.3.5',
        'plotly>=5.14.1',
        'ete3==3.1.2',
        'jellyfish',
<<<<<<< HEAD
        'DendroPy==4.5.2',
        'PyQt5==5.15.9',
        'ray>=2.3.1',
        'deprecated>=1.2.13',
        'psutil==5.9.1',
        'scikit-learn>=1.1.1',
        'scipy>=1.10.1',
        'pyahocorasick>=1.4.4',
=======
        'DendroPy',
        'pyahocorasick',
        'PyQt5',
        'ray',
        'deprecated',
>>>>>>> fcf3bbb741681f8b5cff53167783fe1fd1da4deb

    ],

    entry_points={
        'console_scripts': [
            'cladeomatic=cladeomatic.main:main',
        ],
    },
)