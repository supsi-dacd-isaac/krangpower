PackedOpendssElement reference
''''''''''''''''''''''''''''''
:code:`PackedOpendssElements` are returned by bracket-indicizing a :code:`Krang` with a string containing the name of the desired
element. Aside from electrical components such as :code:`Vsource`, :code:`Line`, etc., buses can also be returned as
:code:`PackedOpendssElement`.
The main objective of :code:`PackedOpendssElement` is to conveniently pack in one object all the methods for the components
:code:`CktElement` interface and the :code:`Vsource` interface; furthermore, in order for the data to be accessed, they have to be
selected via the :code:`Circuit.ActiveElement` method, etc., with one "selector" for each interface. With :code:`PackedOpendssElements`
all the available methods are directly accessible, with no interface swap or the need to explicitly select it with
another method.
In addition to this core functionality, :code:`PackedOpendssElement` enables you to access the underlying object's parameter
values through bracket indexing and also has a number of utility methods that you can see in the class detail docs in
this page.

Class reference
...............

.. autoclass:: krangpower.enhancer._PackedOpendssElement
   :members:
   :special-members:
   :exclude-members: __weakref__