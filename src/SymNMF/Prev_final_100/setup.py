from setuptools import Extension, setup

module = Extension("prev_symnmfmodule", sources=['prev_symnmfmodule.c', 'prev_symnmf.c'], include_dirs=['./'])
setup(name='prev_symnmfmodule',
      version='1.0',
      description='Python wrapper for custom C extension',
      ext_modules=[module])