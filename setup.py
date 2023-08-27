from setuptools import setup, find_packages
from pathlib import Path

# see https://packaging.python.org/guides/single-sourcing-package-version/
version_dict = {}
with open(Path(__file__).parents[0] / "chempotpy/_version.py") as fp:
    exec(fp.read(), version_dict)
version = version_dict["__version__"]
del version_dict

setup(
    name="chempotpy",
    version=version,
    author="Yinan Shu, Zoltan Varga, Donald G. Truhlar",
    description="chempotpy: CHEMical library of POTential energy surfaces in PYthon",
    python_requires=">=3.8",
    packages=find_packages(),
    package_data={
        'chempotpy': ['*/*']
    },
    install_requires=[
        'numpy',
        'gfortran',
    ],
    zip_safe=True,
)
