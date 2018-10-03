#! /usr/bin/env python3

package = 'mageck-ibar'
version = '0.1.1'

import os
import sys
from distutils.command.install import install as DistutilsInstall
from distutils.command.build_py import build_py
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

def readme():
    with open(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'README.md'
            )
    ) as f:
        return f.read()


def compile_rra():
    os.chdir('rra')
    os.system('make')
    rev=os.system('../bin/RRA')
    os.chdir('../')
    return rev


class RRAInstall(DistutilsInstall):
    def run(self):
        # compile RRA
        if(compile_rra()!=0):
            print(
                "CRITICAL: error compiling the RRA source code. Please check your c++ compilation environment.",
                file=sys.stderr
            )
            sys.exit(1)
        DistutilsInstall.run(self)


setup(
    name=package,
    version=version,
    description="The mageck-ibar software is used for analysis of CRISPR-iBAR screening data.",
    long_description=readme(),
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
    ],
    url='',
    author='Zhiheng Liu',
    author_email='zhiheng.liu@pku.edu.cn',
    license='GPL',
    packages=['mibar'],
    install_requires=[
        'numpy', 'scipy', 'pandas'
    ],
    scripts=['bin/mageck-ibar'],
    package_dir={'mibar':'mibar'},
    data_files=[('bin', ['bin/RRA'])],
    cmdclass={'install': RRAInstall, 'build_py': build_py},
    include_package_data=True,
    zip_safe=False
)
