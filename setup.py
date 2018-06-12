from distutils.core import setup
import os

setup(
    name='Robust and sensitive tools for high-throughput data analysis',
    version='1.0',
    packages=['rosely',],
    package_data={'rosely':[os.path.join('rosely', 'HelveticaCY.dfont')]},
    license='Creative Commons',
    long_description=open('README.md').read(),
)
