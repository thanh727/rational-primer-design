from setuptools import setup, find_packages

# Read the README for the long description on PyPI
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="rational-primer-design",
    version="1.0.0",
    author="Thanh Nguyen",
    author_email="nguyenthanh727@gmail.com",
    description="A high-throughput rational design pipeline for TaqMan assays in variable genomes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/rational-primer-design",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.9',
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "primer3-py>=0.6.1",
        "numpy>=1.20.0",
        "tqdm>=4.60.0",
        "python-Levenshtein>=0.12.2"
    ],
    entry_points={
        'console_scripts': [
            'rational-design=rational_design.cli:main',  
        ],
    },
    include_package_data=True,
)