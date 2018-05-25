Curves and dependencies
-----------------------
These are auxiliary object that are used to manage data-points, mainly electrical characteristics and time-dependent
quantities, for sequential simulations (duty-cycle, etc.) of a Circuit. Typical examples are inverter efficiencies,
memorized in Curve (with 'xycurve' parameter), temperature profiles (Curve with 'thsape' parameter) and Load
histories, memorized in loadshapes.

Loadshape
'''''''''
.. autoclass:: krangpower.Loadshape
   :members:
   :undoc-members:
   :inherited-members:

CsvLoadshape
''''''''''''
.. autoclass:: krangpower.CsvLoadshape
   :members:
   :undoc-members:
   :inherited-members:

XYCurve
''''''''''''
.. autoclass:: krangpower.XYCurve
   :members:
   :undoc-members:
   :inherited-members: