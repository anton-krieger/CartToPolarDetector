import io
import os
import sys
from setuptools import setup

# Package meta-data.
NAME = 'CartToPolarDetector'
DESCRIPTION = 'A package to convert a cartesian detector grid into a polar detector grid taking into account the various possible cross section areas of cartesian and polar pixels.'
URL = 'https://github.com/anton-krieger/CartToPolarDetector'
EMAIL = 'akrieger@astrophysik.uni-kiel.de'
AUTHOR = 'Anton Krieger'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.0.1'

# find current path
here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION



setup(
	name=NAME,
	version=VERSION,
	description=DESCRIPTION,
    author=AUTHOR,
    author_email=EMAIL,
	url=URL,
	py_modules=[NAME],
	python_requires=REQUIRES_PYTHON,
	install_requires=[
		"numpy",
	],
	package_dir={'': 'src'},
	license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        "Operating System :: OS Independent",
	],
	
)
