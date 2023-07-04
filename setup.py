#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages
from src.version import __version__
author = 'James Robertson'

classifiers = """
Development Status :: 3 - Beta
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
    name='profile_dists',
    include_package_data=True,
    version=__version__,
    python_requires='>=3.8.2,<4',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/profile_dists',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@phac-aspc.gc.ca',
    description=(
        'Profile Dists: Rapid calcualtion of allele profile distances and distance base querying'),
    keywords='cgMLST, wgMLST, outbreak, surveillance, clustering, distance matrix',
    classifiers=classifiers,
    package_dir={'profile_dists': 'src'},
    package_data={
        "": ["*.txt"],
    },

    install_requires=[
        'pyarrow==12.0.0',
        'fastparquet==2023.4.0',
        'numba==0.57.1',
        'numpy==1.24.4',
        'tables==3.8.0',
        'six>=1.16.0',
        'pandas==2.0.2 ',

    ],

    entry_points={
        'console_scripts': [
            'profile_dists=src.dist:main',
        ],
    },
)