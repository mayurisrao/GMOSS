from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyGMOSS",
    version="0.0.1",
    description="GMOSS (Global MOdel for the radio Sky Spectrum)",
    py_modules=["pyGMOSS"],
    package_dir={"": "pyGMOSS"},
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    author="",
    author_email="",
)
