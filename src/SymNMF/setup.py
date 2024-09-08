from setuptools import Extension, setup

module = Extension("symnmfmodule", sources=['symnmfmodule.c', 'symnmf.c'], include_dirs=['./'])
setup(name='symnmfmodule',
     version='1.0',
     description='Python wrapper for custom C extension for SymNMF',
     ext_modules=[module])

