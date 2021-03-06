.. krangpower documentation master file, created by
   sphinx-quickstart on Fri May 18 12:01:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Krangpower
==========

Krangpower is a Python package for electrical network simulations based on OpenDSS_ by EPRI and a Python api for its direct library incarnation, `OpenDSSDirect.py`_, by Dheepak Krishnamurthy.


Design goals
''''''''''''

- Providing an even easier and more intuitive frontend;
- Introducing **measurement units** through the package pint_, allowing the user to worry less about errors, to demand the burden of conversion and to correctly interpret the results without recurring to the OpenDSS docs;
- Returning results in interoperable data structures and containers;
- Enabling advanced analysis modes, for example:
   - Exporting the network topology in a networkx_ **graph**
   - Solving duty cycles that involve **smart components** that need to update themselves after every n steps
- Providing a **I/O facility** based on json files, allowing easier bookkeping and search and procedural generation of circuits with custom tools


.. _OpenDSS: http://smartgrid.epri.com/SimulationTool.aspx
.. _networkx: https://networkx.github.io/
.. _pint: https://pint.readthedocs.io/
.. _`OpenDSSDirect.py`: https://nrel.github.io/OpenDSSDirect.py/index.html


Installation
''''''''''''

Run as appropriate for your system:

.. code::

   pip install krangpower


.. code::

   pip3 install krangpower

.. important::

   Krangpower is a python 3.x package; currently, compatibility with python 2.x is not sought after.

Alternatively, if you wish to contribute to the codebase too, you can clone and use the `Github repo`_. As of june 2018, krangpower is in intense development and the Github repo is typically a few days/weeks ahead of the pypi index.


.. _`Github repo`: https://github.com/supsi-dacd-isaac/krangpower


Contents
========

.. toctree::
   :maxdepth: 2

   usage
   example_smart
   reference_components
   io
   iockt
   iodss
   oenh_ref
   graph
   defaults_help
   config
   log
   faq
   about

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

