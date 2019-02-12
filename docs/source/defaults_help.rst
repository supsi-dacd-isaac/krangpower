Data files
==========

Krangpower depends on a series of data/configuration files, contained in the :code:`/defaults` folder, in order to offer its functionality.

.. warning::
    These files should only be touched if you know what you're doing, and modifying them improperly could result in unpredictable and wrong behavior.

.. important::

    As of v0.2.2, these files are not yet complete. In order to contribute to the perfectioning of these files, a detailed knowledge of OpenDSS and meticulous study of its documentation is required.
    Completing these configuration files is not the current #1 priority, but it will be in the future, when all the important functionality is completely ironed out, and certainly before any development-state changing release.

    In the meantime, any help on completing (and double-checking!!) these files is appreciated!

In order to aid the knowledgeable user to comprehend what these files do exactly, here follows a brief explanation of
each one.

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

dangerous_stacks.txt
....................

It lists what are the OpenDSS calls that are unsafe to make when a circuit is not yet solved. These calls have an extremely high
risk of retrieving non-physical results or, worse, cause a segmentation fault.
Calls listed here are met with an UnsafeCallError raise instead of being performed when the circuit is unsolved; or, if the "force_unsafe_calls" setting is set to True,
a WARNING is issued and the call is performed anyways.

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


error_strings.json
..................
This file contains those substrings that characterise messages of error coming from the OpenDSS text interface. Since
in these cases OpenDSS does not raise any exception on its own, it's necessary to check the returned strings to see if
they signal an error so that :code:`krangpower.enhancer` can raise an :code:`OpenDSSTextError`, thus avoiding dangerous
and difficult-to-debug silent errors.

.. important::

    If you know of error strings returned by the text interface that are not caught in the cases listed in the file,
    please drop a couple of lines to the developers!

error_strings.json contains two lists:

   - **beginning** is the list of substrings that are found at the beginning of the error message;
   - **middle** is the list of substrings that are found somewhere in the middle of the error message.

mandatory_unit_dec.json
.......................

In vanilla OpenDSS, some element definitions allow you to specify in what units you are going to provide the numerical data (e.g., Lines).
Since in krangpower, physical quantities are provided to the element constructors complete with a pint measurement unit, this specification
becomes superfluous. This file contains the unit declarations chosen by krangpower for its internal workings;
and krangpower then converts what the user provides to these units when generating an OpenDSS instruction.

bypassable_cmd_strings.txt
..........................

This file contains regular expressions that krangpower tries to match to any command sent through :code:`Krang.command()`.
If the command matches any of these regular expressions, the cached :code:`Krang.graph()` (whose calculation can be very expensive)
is mantained.

In other words, this file contains regular expressions that identify those commands that have no influence on the graph and
should not trigger its re-calculation.

splash.txt
..........

Just a splash screen with the logo.