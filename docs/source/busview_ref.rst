BusView reference
'''''''''''''''''
:code:`BusView` are returned by bracket-indicizing a :code:`Krang` with a tuple containing the id of the desired buses.
The main functionality of :code:`BusView` is to enable the possibility of adding to a :code:`Krang` circuit components tied to a
particular set of buses. The operator used to carry out this operation is :code:`lshift`, or :code:`<<`, in the same fashion as direct
addition to a :code:`Krang`.
Aside from this, :code:`BusViews` can also return a list of the components already tied to their set of buses, a list of which
is available in the property content.


Class reference
...............

.. autoclass:: krangpower._krangsuit._BusView
   :members:
   :special-members:
   :exclude-members: __weakref__