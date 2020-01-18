#!/usr/bin/env python
import setuptools
from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils import log

packages = find_packages(exclude=('tests', 'doc'))
provides = ['taurex_ggchem', ]


requires = []


install_requires = ['numpy',
                    'taurex']


def build_fortran(package_name, fortran_sources, extra_compile_args=['-O3','-fdefault-real-8', '-fdefault-double-8','-fbounds-check' ,'-fbacktrace']):
    ext = Extension(name=f'taurex_ggchem.external.{package_name}',
                    extra_compile_args=extra_compile_args,
                    sources=fortran_sources)

    return ext

def build_ggchem():

    sources = [         
                        

                        #'src16/equil_cond.f',
                        
                        'src16/is_nan.F',
                        'src16/nasa_polynomial.f',
                        'src16/ggchem.f',
                        'src16/gauss16.f',
                        'src16/gauss8.f',
                        'src16/gauss_nm.f',
                        'src16/massfrac.f',
                        'src16/linpack_q.f90',
                        'src16/slatec_routines.f',
                        'src16/smchem16.f',
                        'src16/smchem8.f',
                        'src16/supersat.f',
                        'src16/stindex.f',
                        'src16/upper.f',
                        'src16/nucleation.f',
                        'taurex_ggchem/glue/datamod.f90',
                        'taurex_ggchem/glue/taurex_glue.f90',
                        'taurex_ggchem/glue/taurex_glue.pyf',
                        ]

    return build_fortran('fort_ggchem', sources)





entry_points = {'taurex.plugins': 'ggchem = taurex_ggchem'}


extensions = [build_ggchem(), ]
setup(name='taurex_ggchem',
      author='Ahmed Faris Al-Refaie',
      author_email='ahmed.al-refaie.12@ucl.ac.uk',
      license="BSD",
      description='TauREx 3 retrieval framework',
      packages=packages,
      include_package_data=True,
      entry_points=entry_points,
      provides=provides,
      requires=requires,
      install_requires=install_requires,
      ext_modules=extensions
      )