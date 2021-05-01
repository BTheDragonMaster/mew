from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'MEW: mRNA Expression Wizard'
LONG_DESCRIPTION = 'A package to train algorithms that predict mRNA expression level from sequence.'

setup(
    name="mew",
    version=VERSION,
    author="Barbara Terlouw",
    author_email="barbara.terlouw@wur.nl",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=[],
)
