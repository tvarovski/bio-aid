from setuptools import setup, find_packages
import codecs
import os


VERSION = '0.3.2'
DESCRIPTION = 'Genetic Analysis Tools'

with codecs.open("README.md", encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

# Setting up
setup(
    name="bio-aid",
    version=VERSION,
    author="tvarovski (Jerzy Twarowski)",
    author_email="tvarovski1@gmail.com",
    url="https://github.com/tvarovski/bio-aid",
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
