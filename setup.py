from setuptools import setup, find_packages

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(name='subtil',
	version='0.1.0',
	description='Subgroup discovery platform for examination of Bacterial TIme-Lapse microscopy dataset',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://dtai.cs.kuleuven.be/software/python-subtle/',
	author='Antoine Adam',
	author_email='antoine.adam@cs.kuleuven.be',
	license='MIT',
	classifiers=[
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 4 - Beta',

    # Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',

    # Pick your license as you wish (should match "license" above)
     'License :: OSI Approved :: MIT License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
#    'Programming Language :: Python :: 3',
#    'Programming Language :: Python :: 3.2',
#    'Programming Language :: Python :: 3.3',
#    'Programming Language :: Python :: 3.4',
	],
	keywords='time-lapse-microsopy subgroup-discovery microbiology',
	project_urls={
    'Documentation': 'https://dtai.cs.kuleuven.be/software/python-subtle/',
    'Source': 'https://github.com/AntoineAdam/subtil',
    'Tracker': 'https://github.com/AntoineAdam/subtil/issues',
	},
	packages=find_packages(),
	install_requires=['numpy','scipy','matplotlib','scikit-learn'],
	python_requires='>=2.7, <3',
	zip_safe=False,
	#data_files=[("subtil",["subtil/subtildb-structure-sqlite.sql"])]  # for distutils
	package_data={'subtil':["*.sql"]}
	)
