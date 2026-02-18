from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "reweight",
        ["reweight.cpp", "ranlxs.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["-O3"],
    )
]

setup(name="reweight", version="0.1", ext_modules=ext_modules)
