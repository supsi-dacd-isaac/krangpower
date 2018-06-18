Usage
=====

The central class of the krangpower package is the eponymous :code:`Krang`. It is designed to provide easy and intuitive interaction. Let's see the basics of :code:`Krang` usage with a small walkthrough.

.. code-block:: python
   :linenos:

   >>> import krangpower as kp
   >>> um = kp.UM
   >>> src = kp.Vsource(basekv=15000*um.V)
   >>> myKrang = kp.Krang('mighty_krang', src)

- **@2** Krangpower uses the physical quantity management package pint_, so the first thing we did here is to bring in the local namespace krangpower's :code:`UnitRegistry`: the constant :code:`UM`. Think of it as an object whose attributes are the measurement units.
- **@3** Every circuit has to be initialized with a slack voltage source, so :code:`Krang` will need one. In krangpower, voltage sources are represented by :code:`Vsource`. You can initialize :code:`Vsource`, or any other electrical entity, with a set of keyword arguments to override the default values. We will call them "parameters". Here we edited just the parameter 'basekv' to a value of 15000 Volts. More on parameters later.
- **@4** Here we instantiate a :code:`Krang` with a name string and our :code:`Vsource`. Both parameters are optional; if not provided, they will default to :code:`'Unnamed_Krang'` and :code:`kp.Vsource()`. We can now start to add stuff to :code:`myKrang`.

.. code-block:: python
   :linenos:
   :lineno-start: 5

   >>> lc = kp.LineCode_S('lcs1', r0=1*um.ohm/um.km)
   >>> myKrang << lc
   >>> myKrang['sourcebus', 'alpha']
   <BusView('sourcebus', 'alpha')>
   >>> myKrang['sourcebus', 'alpha'] << kp.Line(length=45*um.ft).aka('line1') * lc
   >>> myKrang['line1']
   <PackedOpendssElement: line1>
   >>> myKrang['line1']['r0']
   <Quantity(0.001, 'ohm / meter')>

More irons in the fire.

- **@5** We instantiate a :code:`LineCode_S`. The syntax is the same as the :code:`Vsource`, except that we must specify a name as first parameter. Because they're not made of iron & copper but of information, LineCodes have a name of their own and do not have to be aliased every time they're used (see line 9).
- **@6** Linecodes are meant to be associated to Lines. To make use of :code:`lc`, we first have to inform :code:`myKrang` of its existence. A Linecode is not bound to a bus or couple of buses; so we add it directly. In krangpower, the operator that adds a component to a :code:`Krang` is **left shift**, :code:`<<`.
- **@7** This line demonstrates `bracket indexing Krang with a tuple`_. The object returned is a :code:`BusView`. You can add object that require topological collocation to :code:`BusView`, in the same way of :code:`Krang` themselves.
- **@9** This line has a lot going on.
   - We are instantiating a :code:`Line` specifying its length in feet - we can use any length unit, krangpower will take care of the conversion. In order to be added to the circuit, the line has to be *aliased* with the method :code:`aka`; in this way, it will have a name inside the circuit.
   - The line is then *associated* with :code:`lc`, so it will make use of it in defining its impedance; the *association* operator is multiplication, :code:`*`.
   - The line is added to the :code:`BusView` we saw at point 7.

  Note: The default name of the output bus of the circuit's slack :code:`Vsource` is 'sourcebus'. You can customize it by passing it as a third parameter to :code:`Krang` constructor.
- **@10** Now that 'line1' is inserted into :code:`myKrang`, we can retrieve it by `bracket indexing Krang with a string`_ with its name. The object retrieved is a PackedOpendssElement_.
- **@12** Suffice to say, for now, that we can retrieve a PackedOpendssElement_'s parameters with bracket indexing; here we print its :code:`r0`. As expected, it's a pint_ :code:`Quantity`, and it's coincident with :code:`lc`

.. code-block:: python
   :linenos:
   :lineno-start: 14

   >>> myKrang['alpha']
   <PackedOpendssElement: alpha>
   >>> myKrang['alpha',]
   <BusView('alpha',)>
   >>> myKrang['alpha',] << kp.load(kv=15*um.kV, kW=130*um.kW)

- **@14,@16** Demonstrate again, in the tricky case of buses, bracket indicization with strings and tuples. In the first case, as expected, we get the bus 'alpha' as a :code:`PackedOpendssElement`; in the second case, we get a :code:`BusView` of the single bus 'alpha'.
- **@18** We can add components bound to a single bus to the :code:`BusView`; in this case, for example, we add a :code:`Load`, instantiated with the usual keyword arguments passed to the constructor.


.. toctree::
   :hidden:
   :maxdepth: 1

   krang_ref
   busview_ref
   packed_ref



.. _pint: https://pint.readthedocs.io/
.. _PackedOpendssElement: packed_ref.html
.. _`bracket indexing Krang with a string`: packed_ref.html
.. _`bracket indexing Krang with a tuple`: busview_ref.html