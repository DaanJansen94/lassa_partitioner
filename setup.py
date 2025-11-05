from setuptools import setup, find_packages

setup(
    name="lassa_partitioner",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.80,<=1.90",
    ],
    entry_points={
        'console_scripts': [
            'LassaPartitioner=lassa_partitioner.cli:main',
        ],
    },
    python_requires='>=3.6,<4.0',
    include_package_data=True,
    author="Daan Jansen",
    description="A bioinformatics tool for creating APOBEC3 and non-APOBEC3 partitions from Lassa virus sequence alignments",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/DaanJansen94/lassa_partitioner",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

