[DEV] Data files
================

Krangpower depends on a series of data/configuration files in order to offer its functionality. These files are contained
in the :code:`/defaults` fonder and should only be
touched if you know what you're doing. In order to comprehend what these files do, here follows a brief explanation of
each one.

**As of v0.1.7, these files are not yet complete.** In order to contribute to the perfectioning of these files, a detailed knowledge of opendss and meticulous study of the opendss documentation is required.
Completing these configuration files is not the current #1 priority, but it will be in the future, when all the important functionality is completely ironed out, and certainly before any development-state changing release.

**In the meantime, any help on completing (and checking!!) these files is appreciated!**

association_types.csv
.....................

This csv file contains rows of four elements and its purpose is to decide what happens when the *association* operator
:code:`*` is used to associate one element to the other. For example, the line:

.. code::

    line,linegeometry,geometry,name

means that when one multiplies a :code:`Line` by a :code:`Linegeometry`, the underlying opendss property "geometry" of the Line will be
filled with the Linegeometry's property "name".

In other words, association_types.csv is the file that manages the :code:`*` associations behavior.


default_dssentities.json
........................

This file is nevralgic to :code:`krangpower` component management. It somewhat resembles the lists of components found in
the `JSON i/o files`_. This file defines a "default_<obj>" for each <obj> supported. Inside it, the following dicts can be
found:

**"type"**

It's simply the type of the element. It's there just for practical reasons.

**"properties"**:

They are the default values for the properties, keyed by name.

Implicitly, they define what is the correct data type for that parameter: :code:`string`, :code:`int`, :code:`float`, :code:`list` (:code:`[...]`) or :code:`matrix` (:code:`[[...], [...]]`)

**"units"**

The keys of this dict are a subset of *properties*'.
It contains the measurement unit for the properties internally used by OpenDSS, to which values have to be converted when generating
opendss code and when receiving pint quantities as user input. The measurement units are specified as pint-parsable strings. For matrices with non-homogeneous measurement units,
units in a list can be specified (e.g., linegeometry's x and h)

**"associated"**

It works in the same way as *properties*, but it contains those ones that reflect the association of another object, rather than an intrinsic property (e.g., "geometry" for a Line)

**"topological"**

It works in the same way as *properties*, but it contains those ones that reflect the topological collocation of the object, rather than an intrinsic property (e.g., "bus1" and "bus2" for a Line)


default_settings.json
.....................

This file is somewhat similar to default_dssentities.json. It's used in :code:`Krang` methods :code:`set()` and :code:`get()`. It contains two dicts and a list:

**"values"**

The default value for every opendss option is specified, and, implicitly, its data type.

**"units"**

The keys of this dictionary are a subset of *values*' ones. For the properties listed, the measurement unit for the properties internally used by OpenDSS is indicated

**"contingent"**

A list of those settings that do not influence the identity and physical characterization of a circuit, but are just computational auxiliary parameters, not to be saved in output files.

interf_sel.json
...............

It's used when creating PackedOpendssElement_.

For each object type supported by :code:`krangpower.enhancer`, it lists in **"interfaces"** what are the interfaces in which methods for that kind of
object can be found, and in **selectors** what are the functions to be used to select that object in order to retrieve the correct
data from the interfaces.

interfaces.json
...............

It lists what are the available methods of each submodule of :code:`OpenDSSDirect.py`. It's used in the automatic population
of :code:`krangpower.enhancer` that is performed on load.

**This file is complete**.

measurement_units.json
......................

This file contains, for each :code:`OpenDSSDirect.py` function, what are the units the returned value has, if any. This
information is used when :code:`krangpower.enhancer` wraps a call to :code:`OpenDSSDirect.py` in order to enhance the data
retrieved.

treatments.json
...............

This file contains, for each :code:`OpenDSSDirect.py` function, what chain of functions :code:`krangpower.enhancer` has
to apply to the returned value, if any, in order to enhance the data retrieved. An explanation on what each of the functions
specified does can be found directly in :code:`krangpower.enhancer` code.

.. _PackedOpendssElement: packed_ref.html
.. _`JSON i/o files`: io.html