import setuptools
from pybind11.setup_helpers import Pybind11Extension


ext_modules = [
    Pybind11Extension(
        "count",
        sources=['src/count.cpp']
    )
]


setuptools.setup(
    name="basecount",
    author="Thomas Brier",
    version="1.0.0",
    packages=setuptools.find_packages(),
    ext_modules=ext_modules, # type: ignore
    entry_points = {
        'console_scripts': 'basecount = src.main:run'
    }
)