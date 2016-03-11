from setuptools import setup

from rabifier import __version__, __author__, __title__, __license__


def readme():
    with open('README.rst') as fin:
        return fin.read()

setup(
    name=__title__,
    version=__version__,
    description='A bioinformatic classifier of Rab GTPases',
    long_description=readme(),
    url='https://github.com/evocell/rabifier',
    author=__author__,
    author_email='jarek.surkont@gmail.com',
    packages=['rabifier'],
    license=__license__,
    zip_safe=False,
    keywords=['Rab', 'classifier', 'bioinformatics'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
    ],
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
