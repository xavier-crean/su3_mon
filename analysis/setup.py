from setuptools import setup, Extension
import pybind11
import numpy

ext_modules = [
    Extension(
        "reweight",
        ["reweight.cpp", "ranlxs.cpp"],
        include_dirs=[pybind11.get_include(),
                      numpy.get_include()
                      ],
        language="c++",
        extra_compile_args=["-O3",
                            "-std=c++17",
                            "-D_GNU_SOURCE",
                            ],
        extra_link_args=[
            "-lm"
        ],
    )
]

setup(name="reweight", version="0.1", ext_modules=ext_modules)
