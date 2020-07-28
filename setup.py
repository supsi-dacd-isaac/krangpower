from setuptools import setup
import os
from codecs import open
from docs.source.conf import release

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.md')) as f:
    long_description = f.read()

# with open(os.path.join(here, 'docs', '_version.py'), encoding='utf-8') as f:
#     version = f.read()

# subversioning (usage with pya.bat)
with open(os.path.join(here, 'subversion.txt')) as f:
    subversion = f.read()

# requirements from txt
with open(os.path.join(here, 'requirements.txt')) as f:
    reqs = f.read().splitlines(0)

version = release + subversion

setup(
    name='krangpower',

    version=version,

    description='Distribution System Simulator based on OpenDSS',
    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://krangpower.readthedocs.io',
    download_url='https://github.com/supsi-dacd-isaac/krangpower',

    # Author details
    author='Federico Rosato',
    author_email='federico.rosato@supsi.ch',

    license='MIT',

    packages=['krangpower', 'krangpower.enhancer', 'test', 'krangpower.gv',
              'test.eu_lv', 'test.usage_1', 'test.example_smart'],  # find_packages(),

    install_requires=reqs,

    extras_require={
        "extras": [
        ],
        "dev": [
        ]
    },

    keywords=['OpenDSS', 'pint', 'DSS'],
    py_modules=['krangpower', 'eu_lv_fulltest'],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If us Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'krangpower': [
            'defaults/association_types.csv',
            'defaults/bypassable_cmd_strings.txt',
            'defaults/dangerous_stacks.txt',
            'defaults/default_dssentities.json',
            'defaults/default_settings.json',
            'defaults/DSSHelp.cfg'
            'defaults/error_strings.json',
            'defaults/interf_sel.json',
            'defaults/interfaces.json',
            'defaults/mandatory_unit_dec.json',
            'defaults/measurement_units.json',
            'defaults/splash.txt',
            'defaults/treatments.json',

            'config/krang_config.cfg',

        ]
    },

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',


        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
