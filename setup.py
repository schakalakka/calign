from distutils.core import setup, Extension

# the c extension module
extension_mod = Extension(name="ext_calign", sources=["src/ext_calign.c", "src/calign.c"],
                          extra_link_args=[],
                          extra_compile_args=['-fPIC', '-O3', '-fno-omit-frame-pointer'],
                          )

setup(name="ext_calign", ext_modules=[extension_mod])
