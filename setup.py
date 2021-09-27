#!/usr/bin/env python

from setuptools import setup

setup(name='webnma',
      version='3.4',
      description='Fast normal mode calculation and analyese for protein structures;',
      url='https://github.com/reuter-group/webnma3',
      author='Dandan Xue, Reuter Lab',
      author_email='dandan.xue@uib.no',
      license='',
      packages=['webnma', 'webnma.utils'],
      zip_safe=False,
      ## use conda environment set-up instead of pip
      # install_requires=[
      #     'numpy==1.16.2',
      #     'scipy==1.2.1',
      #     'biopython==1.73',
      #     'matplotlib==3.0.3'
      # ],
      entry_points= {
          'console_scripts': ['webnma=src.cmd_line:main'],
      },
)
