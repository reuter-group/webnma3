#!/usr/bin/env python

from setuptools import setup

setup(name='webnma',
      version='3.4',
      description='Fast normal mode calculation and analyese for protein structures;',
      url='https://github.com/reuter-group/webnma3',
      author='Dandan Xue, Reuter Lab',
      author_email='dandan.xue@uib.no',
      license='GPL-3.0 License',
      package_dir = {'webnma': 'src'},
      packages=['webnma', 'webnma.utils'],
      zip_safe=False,
      ## recommand to use conda environment to install the dependencies
      # install_requires=[
      # ],
      entry_points= {
          'console_scripts': ['webnma=webnma.cmd_line:main'],
      },
)
