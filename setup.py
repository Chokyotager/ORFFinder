#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="orffinder",
    version="1.4",
    author="ChocoParrot",
    author_email="lachocoparrot@gmail.com",
    description="ORFFinder API.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Chokyotager/ORFFinder",
    project_urls={
        "Bug Tracker": "https://github.com/Chokyotager/ORFFinder/issues",
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6"
)
