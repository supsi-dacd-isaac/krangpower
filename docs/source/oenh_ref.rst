The krangpower.enhancer submodule
'''''''''''''''''''''''''''''''''

The submodule :code:`krangpower.enhancer` constitues an `OpenDSSDirect.py`_ wrapper and overhaul, and is designed to expose the exact same interface - the api reference is the same you can find in its docs. It is usable in and by itself as a substitute of the `OpenDSSDirect.py`_ module when interested in the additional functionality.

This means that if:

.. code::

   opendssdirect.dss.<X>

is a valid expression, then:

.. code::

   krangpower.ehnancer.<X>

also is, but the results returned by :code:`krangpower.enhancer` have several advantages:

- Items that OpenDSS returns as simple lists of floats (e.g., :code:`opendssdirect.dss.Circuit.Losses()`, :code:`opendssdirect.dss.Bus.Voltages()`, etc.) are returned as lists of complex numbers, matrices of [nterm x ncond], etc. as is appropriate.
- Structured items such as :code:`opendssdirect.Circuit.SistemY()` are returned as :code:`pandas.DataFrame` for easier manipulation and export
- Items come, where appropriate, as Quantities (from the pint_ package) with the appropriate measurement unit. This enables easy conversions and secures against miscalculations.
- Through the exposed function :code:`txt_command`, the OpenDSS text interface is checked for errors (that are normally just returned as strings without raising anything).
- The exposed a :code:`pack` function that enables easy and intuitive element exploration by returning a PackedOpendssElement_ corresponding to the name passed as argument, either fully qualified or simple.

Additional functions
....................

:code:`krangpower.enhancer` exposes a few more utility functions.

.. automodule:: krangpower.enhancer
   :members: pack, get_all_names, txt_command, log_line


.. _PackedOpendssElement: packed_ref.html
.. _pint: https://pint.readthedocs.io/
.. _`OpenDSSDirect.py`: https://nrel.github.io/OpenDSSDirect.py/index.html