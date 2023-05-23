from setuptools import setup, find_packages
import codecs
import os


VERSION = '0.3.0'
DESCRIPTION = 'Genetic Analysis Tools'
LONG_DESCRIPTION = 'A package that contains tools for the analysis of genetic sequence and sequencing related data'

# Setting up
setup(
    name="bio-aid",
    version=VERSION,
    author="tvarovski (Jerzy Twarowski)",
    author_email="<tvarovski1@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['pandas', 'numpy', 'matplotlib', 'seaborn', 'regex', 'pyensembl', 'natsort'],
    keywords=['python', 'biology', 'bio', 'genetics', 'genomics', 'NGS'],
    classifiers=[
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3.10",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
