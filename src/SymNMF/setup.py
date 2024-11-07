from setuptools import Extension, setup

module = Extension("symnmfmodule", sources=['symnmfmodule.c', 'symnmf.c'], include_dirs=['./', './Prev_final_100'])
setup(name='symnmfmodule',
     version='1.0',
     description='Python wrapper for custom C extension for SymNMF',
     ext_modules=[module])

module2 = Extension("prev_symnmfmodule", sources=['Prev_final_100/prev_symnmfmodule.c', 'Prev_final_100/prev_symnmf.c'], include_dirs=['./'])
setup(name='prev_symnmfmodule',
      version='1.0',
      description='Python wrapper for custom C extension',
      ext_modules=[module2])

