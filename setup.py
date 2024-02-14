#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = 'juncmut',
    version = '0.5.3',
    description='Python programs for the identification of the genomic mutation from RNA-seq data',
    url = 'https://github.com/ncc-gap/juncmut.git',
    author = 'Naoko Iida, Ai Okada and Yuichi Shiraishi',
    author_email = 'genomon.devel@gmail.com',
    license = 'GPLv3',

    classifiers = [
        'Development Status :: 5 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: ',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(),
    package_data={},
    #install_requires = [],

    entry_points = {'console_scripts': ['juncmut = juncmut:main']}

)
