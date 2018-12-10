[![Documentation Status](https://readthedocs.org/projects/krangpower/badge/?version=master)](https://krangpower.readthedocs.io/en/master/?badge=master)
[![Build Status](https://travis-ci.org/supsi-dacd-isaac/krangpower.svg?branch=master)](https://travis-ci.org/supsi-dacd-isaac/krangpower)
[![codecov](https://codecov.io/gh/supsi-dacd-isaac/krangpower/branch/master/graph/badge.svg)](https://codecov.io/gh/supsi-dacd-isaac/krangpower)
[![Latest Version](https://img.shields.io/pypi/v/krangpower.svg)](https://pypi.python.org/pypi/krangpower/)

# Krangpower
Distribution System Simulator based on [OpenDSS](https://sourceforge.net/projects/electricdss/) and [OpenDSSDirect.py](https://nrel.github.io/OpenDSSDirect.py/index.html). Modern Syntax, DataFrames, Pint, Networkx.

# Design goals
* Providing an even easier and more intuitive frontend
* Introducing measurement units through the package [pint](http://pint.readthedocs.io/en/latest/), allowing the user to worry less about errors, to delegate the burden of conversion and to correctly interpret the results without recurring to the OpenDSS docs
* Returning results in interoperable data structures and containers
* Enabling advanced analysis modes, for example:
    * Exporting the network topology in a [networkx](https://networkx.github.io/) graph
    * Solving duty cycles that involve smart components that need to update themselves every n steps
* Providing a I/O facility based on zip packages and JSON files, allowing easier bookkeping and search and procedural generation of circuits with custom tools
* And many more!

# Documentation
The documentation for krangpower is available [here](https://krangpower.readthedocs.io).

# Python version
Krangpower is a python 3.x package; currently, compatibility with python 2.x is not sought after.

# Acknowledgements
The authors would like to thank the Swiss Federal Office of Energy (SFOE) and the Swiss Competence Center for Energy Research - Future Swiss Electrical Infrastructure (SCCER-FURIES), for their financial and technical support to this research work.

