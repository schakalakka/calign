from distutils.core import setup, Extension

# the c extension module
extension_mod = Extension(name="swalign", sources=["src/swalign.c"])

setup(name = "swalign", ext_modules=[extension_mod])