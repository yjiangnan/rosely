from distutils.core import setup
import os

setup(
    name='rosely',
    description='Robust and sensitive tools for high-throughput data analysis',
    version='1.2.2',
    packages=['rosely',],
    package_dir={'rosely': 'rosely'},
    package_data={'rosely':['HelveticaCY.dfont']},
    license='Creative Commons',
    long_description=open('README.md').read(),
)
