import os
from setuptools import setup

name = "velocity-estimation"

with open("README.md") as f:
    long_description = f.read()

here = os.path.abspath(os.path.dirname(__file__))

setup(
    name=name,
    description="Python scripts for velocity estimation in turbulent flows",
    author="Juan Manuel Losada",
    author_email="juan.m.losada@uit.no",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/uit-cosmo/velocity-estimation",
    download_url="https://github.com/uit-cosmo/velocity-estimation/archive/refs/tags/1.1.tar.gz",
    license="MiT",
    version="1.0.0",
    packages=["velocity_estimation"],
    python_requires=">=3.0",
    install_requires=[
        "blobmodel>=1.0.0",
        "numpy>=1.23.3",
        "scipy>=1.11.3",
    ],
    classifiers=[
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    zip_safe=False,
)
