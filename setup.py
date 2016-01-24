from distutils.core import setup, Extension

# the c extension module
extension_mod = Extension(name="ext_swalign", sources=["src/ext_swalign.c", "src/swalign.c"],
                          extra_link_args=[],
                          extra_compile_args=['-fPIC', '-O3', '-fno-omit-frame-pointer'],
                          )

setup(name="ext_swalign", ext_modules=[extension_mod])
