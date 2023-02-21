"""
Setup Model
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="AdenineFootprinter",
    version="1.0.0",
    description="Find Nucleosome Footprint Using methyladenosine in genome",  # Optional
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zhuweix/AdenineFootprinter",
    author="Zhuwei Xu",
    author_email="wszwei@gmail.com",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="nucleosome footprint, genome accessibility, methyladenosine, methylated adenine, computational biology",  # Optional
    package_dir={"": "src",
                 "example": "example"},
    package_data={"src": ["asset/*.data"],
                  "example": ["figure/*.png"]},
    packages=find_packages(where="src"),
    python_requires=">=3.7, <4",
    install_requires=["numpy",
                      "scipy",
                      'pandas',
                      "statsmodels",
                      'pysam',
                      'biopython'
                      ],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            # command = package.module:function
            'footprinter = adeninefootprinter.footprinter:main',
        ],},
    project_urls={
        "Bug Reports": "https://github.com/zhuweix/AdenineFootprinter/issues",
        "Source": "https://github.com/zhuweix/AdenineFootprinter",
        "Lab": "https://www.nichd.nih.gov/research/atNICHD/Investigators/clark",
    },
)