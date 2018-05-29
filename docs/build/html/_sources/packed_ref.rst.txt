PackedOpendssElement reference
''''''''''''''''''''''''''''''
:code:`PackedOpendssElements` are returned by bracket-indicizing a :code:`Krang` with a string containing the name of the desired
element. Aside from electrical components such as :code:`Vsource`, :code:`Line`, etc., buses can also be returned as
:code:`PackedOpendssElement`.
The main objective of :code:`PackedOpendssElement` is to conveniently pack in one object all the methods for the components
that are present in the various opendss interfaces (for example, for :code:`Vsource`, from :code:`opendssdirect.CktElement` and :code:`opendssdirect.Vsources`); Furthermore, with the classic opendss interface, before obtaining the correct data, the component has to be
selected via other interfaces, (again, for :code:`Vsource`, :code:`opendssdirect.Circuit.ActiveElement` and :code:`opendssdirect.Vsources.Name`).
With :code:`PackedOpendssElements`
all the available methods are directly accessible, with no interface swap or the need to explicitly select the component
In addition to this core functionality, :code:`PackedOpendssElement` enables you to access the underlying object's parameter
values through bracket indexing and also has a number of utility methods that you can see in the class detail docs in
this page.

Class reference
...............

.. autoclass:: krangpower.enhancer.OpendssdirectEnhancer._PackedOpendssElement
   :members:
   :special-members:
   :exclude-members: __weakref__