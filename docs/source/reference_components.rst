Components
==========

Components are classes exposed by :code:`krangpower` that identify electrical components, monitors, linecodes...
They typically are instantiated with keyword arguments to set their parameters, possibly further modified, and then
added to a Krang with the :code:`<<` operator after, if necessary, associating them to one another with the :code:`*` operator.

the paramhelp() method
......................

Don't forget that every instance of a component has a :code:`paramhelp()` method, that will print to console useful information about parameters.

reference pages
...............

Here you can find reference for each component available in krangpower. In particular, a table with the component's
parameters is provided.

.. toctree::
   :maxdepth: 1

   electrical
   smart
   above
   linedef
   condef
   curves




