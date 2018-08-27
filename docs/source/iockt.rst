I/O with zip archives
=====================

:code:`Krang` allows to export a package that encloses in one place everything needed to reproduce the circuit.
A bare JSON_, in fact, contains all the topological information, but references such as csv loadprofile files, :code:`DecisionModels`
for the smart components and the likes are not included.
In order to save/load a full pack, the following functions can be used:

Saving
......

.. automethod:: krangpower.Krang.pack_ckt
   :noindex:

Loading
.......

.. automethod:: krangpower.Krang.open_ckt
   :noindex:

The files saved are zip archives that contain, aside from the main JSON file, other archives and pickle-files with all the
necessary information.

md5 hash
........

On save and load, the Krang is hashed with the :code:`fingerprint` method in order to verify correctness. The hash
can be seen in a .md5 file present in the ckt archive.


.. _JSON: io.html