Krang reference
'''''''''''''''
The :code:`Krang` is the main class of krangpower. It has facilities to drive the OpenDSS engine, retrieve and query
elements and options, saving and loading circuit data as JSON, representing the circuit as graph, perform repeated
solutions and reporting the results as :code:`DataFrames`. The basic workflow pertaining the :code:`Krang` is illustrated
in the `Usage page`_.

It's worth spending a word about :code:`Krang.brain`. This is an instance of :code:`OpendssdirectEnhancer` and is designed
to expose the exact same interface as the module `OpenDSSDirect.py`_. This means that you can use :code:`Krang.brain`
exactly as you would use this module, with two advantages:

- Results are returned in enhanced data containers as suitable, rather than as raw strings (integers will be of type :code:`int`, arrays will be of type :code:`numpy.ndarray`, pyhysical quantities will have a pint_ unit measure);
- It's bracket-indicizable with a string and returns a :code:`PackedOpendssElement`.

More info on :code:`OpendssdirectEnhancer` is available at the `OpendssdirectEnhancer reference page`_.

Class reference
...............
.. autoclass:: krangpower.Krang
   :members:
   :special-members:
   :exclude-members: __weakref__


.. _`Usage page`: usage.html
.. _`OpenDSSDirect.py`: https://nrel.github.io/OpenDSSDirect.py/index.html
.. _pint: https://pint.readthedocs.io/
.. _PackedOpendssElement: packed_ref.html
.. _`OpendssdirectEnhancer reference page`: oenh_ref.html