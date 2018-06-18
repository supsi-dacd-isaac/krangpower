I/O with JSON files
'''''''''''''''''''
Krang has facilities to store all the relevant information they contain in structured JSON files. The essential structure of such file is this:

.. code::

    krangpower JSON v0.1.6

    root
    ├──"ckt_name": <string>
    ├──"elements":
    │   ├──"type1.name1":
    │   │    ├──"type": <string>
    │   │    ├──"name": <string>
    │   │    ├──"properties": {<prop_name>: <value>}
    │   │    ├──("units"): {<prop_name>: <unit_string>}
    │   │    ├──("depends"): {<prop_name>: <name_string>}
    │   │    ├──("topological"): [<bus_name_string>]
    │   │    └──("hash"): <string> (*)
    │   └──"type2.name2":
    │       '.....
    ├──"settings":
    │   ├──"values": {<setting_name>: <value>}
    │   └──"units": {<setting_name>: <unit_string>}
    └──"buscoords": {<bus_name_string>: [x y]}

    # NOTES:
    #  keys in parentheses are not necessarily present
    #  (*) "hash" is used only for comparisons; when editing a json externally, it can be omitted

Full examples are available under the package's :code:`/test` folder. The methods here illustrated are documented in the `Krang reference page`_ too.

Saving
......

:code:`Krang` can output such a json file through the :code:`Krang.save_json()` method.


Loading
.......

.. autofunction:: krangpower.from_json

A function :code:`from_json(path)` is exposed directly by krangpower. It returns a :code:`Krang` in the state syntesised by the JSON.

Note: inside the box, :code:`from_json(path)` instantiates an empty :code:`Krang`, adds all the components in the right order of dependency and finally returns it.


.. _`Krang reference page`: krang_ref.html