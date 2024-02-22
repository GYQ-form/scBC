from setuptools import setup, find_packages

__version__ = "0.3.2"

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name="scBC",
    version=__version__,
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    install_requires=[
        'scvi-tools>=0.19.0',
        'biomart',
        "requests",
        "scipy"
        ],
    author="Yuqiao Gong",
    author_email="gyq123@sjtu.edu.cn",
    keywords=["single cell transcriptomics", "Biclustering", "bioinformatics", "variational inference(VI)"],
    description="a single-cell transcriptome Bayesian biClustering framework",
    license="MIT",
    url='https://github.com/GYQ-form/scBC',
    long_description_content_type='text/markdown',
    long_description=long_description
)