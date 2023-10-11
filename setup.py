from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="pyGMOSS",
    version="0.0.1",
    description="A Radio sky spectrum generator",
    package_dir={"": "pyGMOSS"},
    packages=find_packages(where="pyGMOSS"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mohithpa/GMOSS/tree/pyGMOSS",
    author="Dr. Mayuri S. Rao, Mohith P A",
    author_email="mohithpa2@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    install_requires=["pandas", "scipy", "numpy", "healpy"],
    extras_require={
        "dev": ["twine>=4.0.2"],
    },
    python_requires=">=3.7",
)