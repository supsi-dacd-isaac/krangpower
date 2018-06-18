Krang reference
'''''''''''''''
The :code:`Krang` is the main class of krangpower. It has facilities to drive the OpenDSS engine, retrieve and query
elements and options, saving and loading circuit data as JSON, representing the circuit as graph, perform repeated
solutions and reporting the results as :code:`DataFrames`. The basic workflow pertaining the :code:`Krang` is illustrated
in the `Usage page`_.

It's worth spending a word about the attribute :code:`Krang.brain`. It points to the :code:`krang.enhancer` submodule and is designed
to expose the exact same interface as the module `OpenDSSDirect.py`_. This means that you can use :code:`Krang.brain`
exactly as you would use this module, with two advantages:

- Results are returned in enhanced data containers as suitable, rather than as raw strings (integers will be of type :code:`int`, arrays will be of type :code:`numpy.ndarray`, pyhysical quantities will have a pint_ unit measure);
- It's bracket-indicizable with a string and returns a :code:`PackedOpendssElement`.

**More info** on :code:`krang.enhancer` is available at the `krang.enhancer reference page`_.

.. IMPORTANT::
   :code:`Krang` is a singleton_, meaning that you cannot have more than one instance around at any given time - an Exception will be raised when trying to instantiate more than one.
   This is because OpenDSS supports only one circuit at any given time. If, in the future, the capability to have more than one circuit will be built in Opendss, this limitation will be removed.

.. NOTE::
   :code:`Krang` can be deleted with the :code:`del` keyword. After doing it, it will be possible (in principle) to instantiate a new one; beware, though, that if you made *strong* references to it
   (for example, assigning it to variables, to lists, to object attributes) the new instantiation will fail even after the :code:`del` command. In other words, a new instantiation will be possible after you del the last
   *strong* reference that was made to the Krang.

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
.. _`krang.enhancer reference page`: oenh_ref.html
.. _singleton: https://en.wikipedia.org/wiki/Singleton_pattern