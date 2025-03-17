from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="plasann",
    version="1.0.7",
    author="Habibul Islam",
    author_email="hislam2@ur.rochester.edu",
    description="A plasmid annotation tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ajlopatkin/PlasAnn",
    project_urls={
        "Bug Tracker": "https://github.com/ajlopatkin/plasann/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=[
        "biopython>=1.78",
        "pandas>=1.0.0",
        "matplotlib>=3.0.0",
        "pycirclize>=0.1.0",
        "gdown>=4.0.0",
    ],
    entry_points={
        "console_scripts": [
            "PlasAnn=Scripts.annotate_plasmid:cli",
        ],
    },
)