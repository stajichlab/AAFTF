import os
from setuptools import setup

from AAFTF.__version__ import __version__
version = __version__

long_description = """
``AAFTF`` Automatic Assembly For The Fungi is a toolkit for automated
          genome assembly, cleanup, mitochondrial genome asm, and polishing'
"""

HERE = os.path.dirname(__file__)

with open(os.path.join(HERE, "requirements.txt"), "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
    name="AAFTF",
    version=version,
    install_requires=install_requires,
    requires=['python (>=3.6.0)'],
    packages=['AAFTF',
              'scripts'],
    author="Jason Stajich, Jonathan Palmer",
    description='Automated genome assembly, cleanup, and polishing',
    long_description=long_description,
    url="http://github.com/stajichlab/AAFTF",
    package_dir={'AAFTF': "AAFTF"},
    package_data={'AAFTF': ['test']},
    zip_safe=False,
    include_package_data=True,
    # scripts=['AAFTF/scripts/AAFTF'],
    entry_points={
        'console_scripts': [
            'AAFTF=AAFTF.AAFTF_main:main',
        ],
    },
    author_email="jasonstajich.phd@gmail.com",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
    )
