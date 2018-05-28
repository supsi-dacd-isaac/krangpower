BusView reference
'''''''''''''''''
BusView are returned by bracket-indicizing a Krang with a tuple containing the id of the desired buses.
The main functionality of BusView is to enable the possibility of adding to a Krang circuit components tied to a
particular set of buses. The operator used to carry out this operation is lshift, or <<, in the same fashion as direct
addition to a Krang.
Aside from this, BusViews can also return a list of the components already tied to their set of buses, a list of which
is available in the property content.

.. autoclass:: krangpower.krangsuit._BusView
   :members:
   :special-members:
