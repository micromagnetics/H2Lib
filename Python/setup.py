from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

ext1 = Extension("lindholm",
                 sources = ["lindholm.pyx", "liblindholm.cxx"],
                 include_dirs = ["../Library"],
                 library_dirs = ["../build/Library"],
                 libraries = ["m", "H2Lib"],
                 extra_compile_args = ["-std=gnu++11", "-O3"],
#                 extra_link_args=['-lgomp'],
                 language="c++"
)

setup(name='lindholm',
      version='0.1',
      description='H2Matrix using Analytic Doublelayer Potential proposed by Lindholm',
      author='Florian Bruckner',
      author_email='florian.bruckner@tuwien.ac.at',
      ext_modules = cythonize(ext1)
)
