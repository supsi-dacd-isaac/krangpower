.. krangpower documentation master file, created by
   sphinx-quickstart on Fri May 18 12:01:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Krangpower
==========

Krangpower is a Python package for electrical network simulations based on OpenDSS_ by EPRI and a Python api for its direct library incarnation, `OpenDSSDirect.py`_, by Dheepak Krishnamurthy.


Design goals:
'''''''''''''

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


Contents tree
=============

.. toctree::
   :maxdepth: 2

   usage
   reference_components
   io
   faq

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

