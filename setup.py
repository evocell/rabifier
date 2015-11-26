from setuptools import setup

from rabifier import __version__, __author__, __title__, __license__


def readme():
    with open('README.md') as fin:
        return fin.read()

setup(
    name=__title__,
    version=__version__,
    long_description=readme(),
    url='https://github.com/evocell/rabifier',
    author=__author__,
    author_email='jarek.surkont@gmail.com',
    packages=['rabifier'],
    license=__license__,
    zip_safe=False,
    install_requires=[
        'biopython',
        'numpy',
        'scipy'
    ],
    setup_requires=[
        'numpy'
    ],
    extras_require={
        'plotting': ['matplotlib']
    },
    include_package_data=True,
    scripts=[
        'bin/rabifier',
        'bin/rabifier-mkdb'
    ]
)
