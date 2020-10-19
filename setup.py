#!/usr/bin/env python

from setuptools import setup, find_packages
 
    
setup(
    name = 'juncmut',
    version = '0.0.2a',
    description='Python programs for the identification of the genomic mutation from RNA-seq splicing junction data',
    url = 'https://github.com/ni6o6/',
    author = 'Naoko Iida',
    author_email = 'iida.nao08@gamil.com',
    license = '',

    classifiers = [
        'Development Status :: 1 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(exclude = ['tests']),
    
    install_requires = [],

    entry_points = {'console_scripts': ['juncmut = juncmut:main']}

)
