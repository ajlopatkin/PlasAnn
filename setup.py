from setuptools import setup, find_packages

setup(
    name='PlasAnn',
    version='1.0.1',
    author='Habibul Islam',
    author_email='hislam2@ur.rochester.edu',
    description='A tool for plasmid annotation and visualization',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas',
        'matplotlib',
        'pycirclize',
        'biopython',
        'gdown',
    ],
    entry_points={
        'console_scripts': [
            'PlasAnn=plasann.annotate_plasmid:main',
        ],
    },
)