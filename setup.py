#!/usr/bin/env python

from setuptools import setup, find_packages

kwargs = {
    'name': 'openmc_uq',
    'version': '0.0.1',
    'packages': find_packages(include=['openmc_uq']),
    'author': 'Ander Gray',
    'author_email': 'ander.gray@ukaea.uk',
    'license' : 'MIT',
    'download_url': 'https://github.com/AnderGray/openmc_UQ',
    'url': 'https://github.com/AnderGray/openmc_UQ',
}

setup(**kwargs)
