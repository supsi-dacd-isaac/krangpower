PackedOpendssElement reference
''''''''''''''''''''''''''''''
PackedOpendssElements are returned by bracket-indicizing a Krang with a string containing the name of the desired
element. Aside from electrical components such as Vsource, Line, etc., buses can also be returned as
PackedOpendssElement.
The main objective of PackedOpendssElement is to conveniently pack in one object all the methods for the components
that are available in the OpenDSS direct interface in various places. Vsources, for example, have methods accessible via both the
CktElement interface and the Vsource interface; furthermore, in order for the data to be accessed, they have to be
selected via the Circuit.ActiveElement method, etc., with one "selector" for each interface. With PackedOpendssElements
all the available methods are directly accessible, with no interface swap or the need to explicitly select it with
another method.
In addition to this core functionality, PackedOpendssElement enables you to access the underlying object's parameter
 values through bracket indexing and also has a number of utility methods that you can see in the class detail docs in
 this page.

.. autoclass:: krangpower.enhancer._PackedOpendssElement
   :members:
   :special-members:
