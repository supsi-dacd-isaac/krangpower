I/O with JSON files
'''''''''''''''''''
Krang has facilities to store all the relevant information they contain in structured JSON files. To save a JSON, the
method used, visible in the `Krang reference page`_, is save_json. The essential structure of such file is this:

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
    │   │    └──("topological"): [<bus_name_string>]
    │   └──"type2.name2":
    │       '.....
    ├──"settings":
    │   ├──"values": {<setting_name>: <value>}
    │   └──"units": {<setting_name>: <unit_string>}
    └──"buscoords": {<bus_name_string>: [x y]}





.. _`Krang reference page`: krang_ref.html