import setuptools
from pybind11.setup_helpers import Pybind11Extension, build_ext


ext_modules = [
    Pybind11Extension(
        "count",
        sources=['basecount/count.cpp']
    )
]


setuptools.setup(
    name="basecount",
    author="Thomas Brier",
    version="1.2.0",
    packages=setuptools.find_packages(),
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules, # type: ignore
    entry_points = {
        'console_scripts': 'basecount = basecount.main:run'
    }
)
