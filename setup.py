import setuptools
from pybind11.setup_helpers import Pybind11Extension, build_ext


ext_modules = [
    Pybind11Extension(
        "count",
        sources=['basecount/count.cpp']
    )
]

exec(open('basecount/version.py').read())

setuptools.setup(
    name="basecount",
    author="Thomas Brier",
    version=__version__, # type: ignore
    packages=setuptools.find_packages(),
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules, # type: ignore
    entry_points = {
        'console_scripts': 'basecount = basecount.main:run'
    }
)
