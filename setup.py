#!/usr/bin/env python

from setuptools import setup, find_packages
from juncmut import __version__

setup(
    name = 'juncmut',
    version = __version__,
    description='Python programs for the identification of the genomic mutation from RNA-seq data',
    url = 'https://github.com/ncc-gap/juncmut.git',
    author = 'Naoko Iida, Ai Okada and Yuichi Shiraishi',
    author_email = 'yuishira@ncc.go.jp, aokada@ncc.go.jp',
    license = 'GPLv3',

    classifiers = [
        'Development Status :: 5 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(exclude = ['resource']),
    package_data={'annot_utils': ['data/*/*']},
    install_requires = [],

    entry_points = {'console_scripts': ['juncmut = juncmut:main', 'annot_utils = annot_utils:main']}
)
