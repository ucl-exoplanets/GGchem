project(
  'TaurexGGChem',
  'c',
  version: '1.1.1',
  license: 'BSD-3',
  meson_version: '>= 1.1.0',
  default_options: [
    'c_std=c99',
    'fortran_std=legacy',
    'optimization=2',
  ],
)

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()

# Emit a warning for 32-bit Python installs on Windows; users are getting
# unexpected from-source builds there because we no longer provide wheels.
is_windows = host_machine.system() == 'windows'
if is_windows and py3.has_variable('EXT_SUFFIX')
  ext_suffix = py3.get_variable('EXT_SUFFIX')
  if ext_suffix.contains('win32')
    warning('You are building from source on a 32-bit Python install. SciPy does not provide 32-bit wheels; install 64-bit Python if you are having issues!')
  endif
endif

cc = meson.get_compiler('c')
# generator() doesn't accept compilers, only found programs - cast it.

# Check compiler is recent enough (see "Toolchain Roadmap" for details)
if cc.get_id() == 'gcc'
  if not cc.version().version_compare('>=8.0')
    error('TaurexGGchem requires GCC >= 8.0')
  endif
elif cc.get_id() == 'msvc'
  if not cc.version().version_compare('>=19.20')
    error('TaurexGGchem requires at least vc142 (default with Visual Studio 2019) ' + \
          'when building with MSVC')
  endif
endif

_global_c_args = cc.get_supported_arguments(
  '-Wno-unused-but-set-variable',
  '-Wno-unused-function',
  '-Wno-conversion',
  '-Wno-misleading-indentation',
)
add_project_arguments(_global_c_args, language : 'c')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

if host_machine.system() == 'os400'
  # IBM i system, needed to avoid build errors - see gh-17193
  add_project_arguments('-D__STDC_FORMAT_MACROS', language : 'cpp')
  add_project_link_arguments('-Wl,-bnotextro', language : 'c')
  add_project_link_arguments('-Wl,-bnotextro', language : 'cpp')
  add_project_link_arguments('-Wl,-bnotextro', language : 'fortran')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

if ff.has_argument('-fdefault-real-8')
  add_project_arguments('-fdefault-real-8', language: 'fortran')
endif
if ff.has_argument('-fdefault-double-8')
  add_project_arguments('-fdefault-double-8', language: 'fortran')
endif

if ff.has_argument('-r8')
  add_project_arguments('-r8', language: 'fortran')
endif


# Intel compilers default to fast-math, so disable it if we detect Intel
_intel_cflags = []
_intel_fflags = []
if cc.get_id() == 'intel'
  _intel_cflags += cc.get_supported_arguments('-fp-model=strict')
elif cc.get_id() == 'intel-cl'
  _intel_cflags += cc.get_supported_arguments('/fp:strict')
endif
if ff.get_id() == 'intel'
  _intel_fflags = ff.get_supported_arguments('-fp-model=strict')
  minus0_arg = ['-assume', 'minus0']
  if ff.has_multi_arguments(minus0_arg)
    _intel_fflags += minus0_arg
  endif
elif ff.get_id() == 'intel-cl'
  # Intel Fortran on Windows does things differently, so deal with that
  # (also specify dynamic linking and the right name mangling)
  _intel_fflags = ff.get_supported_arguments(
    '/fp:strict', '/MD', '/names:lowercase', '/assume:underscore',
    '/assume:minus0'
  )
endif
add_project_arguments(_intel_cflags, language: ['c', 'cpp'])
add_project_arguments(_intel_fflags, language: 'fortran')


subdir('taurex_ggchem/')