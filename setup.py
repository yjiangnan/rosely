from distutils.core import setup
import os

setup(
    name='Robust and sensitive tools for high-throughput data analysis',
    version='1.2.1',
    packages=['rosely',],
    package_dir={'rosely': 'rosely'},
    package_data={'rosely':['HelveticaCY.dfont']},
    license='Creative Commons',
    long_description=open('README.md').read(),
)
