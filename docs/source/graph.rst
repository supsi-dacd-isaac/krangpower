Circuit Graph and Views
=======================

krang.graph
...........

:code:`Krang` has a graph attribute that returns a networkx_ :code:`Graph`. The nodes of this graph are constituted by the buses; the
edges are the elements that connect these buses, such as :code:`Line`, :code:`Transformer`, :code:`Reactor`...

The nodes have:

- a property :code:`'bus'`, which is the :code:`PackedOpendssElement` corresponding to that bus;
- a property :code:`'el'`, which is a list of the :code:`PackedOpendssElement` that are bound to that bus;

The edges have:

- a property :code:`'el'`, which is a list of the :code:`PackedOpendssElement` that are bound to that couple of buses.

Since it contains a plethora of :code:`PackedOpendssElement`, the graph has the potentiality for retrieving all sorts of useful
information and is amenable to data query and reorganization with :code:`networkx`.

GraphView
.........

To aid in working with the graph, krangpower exposes the template class :code:`GraphView`. The :code:`GraphView` is initialized with a
:code:`Krang` instance and two arguments, :code:`busfun` and :code:`edgefun`, that are functions meant to be applied to every :code:`'el'[0]` of the
edges and to every :code:`'bus'` in the nodes (these function can also be :code:`None`).

Initialized in this way, the :code:`GraphView` is then bracket-indicizable with bus names
(for nodes) or couples of bus names (for edges) and retrieves the result of the application of the function on that bus
or edge.

.. autoclass:: krangpower.GraphView
   :members:
   :special-members:

Another method of usage is to inherit :code:`GraphView` and override its __init__ method in order to pack busfun and
edgefun inside it. This is done, for example, with the built-in :code:`VoltageView`, that can be found in the submodule :code:`gv` and initialized with just
a :code:`Krang` instance and returns the property :code:`Voltages()` of both :code:`'bus'` at buses and :code:`'el'[0]` at edges.
More specialized built-in :code:`GraphView` s can and will be added in the future to :code:`krangpower.gv`.

Example:

.. code::

    >>> import krangpower as kp
    >>> from numpy import round
    >>> 
    >>> um = kp.UM
    >>> src = kp.Vsource(basekv = 10.0 * um.kV)
    >>> twc = kp.Krang('twc', src)
    >>> twc['sourcebus', 'a'] << kp.Line(units= 'm', length=20.0 * um.m).aka('l1')
    >>> trf = kp.Transformer(windings=3,
    >>>                     conns=['delta', 'wye', 'wye'],
    >>>                     kvs=[10.0, 2.0, 3.0] * um.kV,
    >>>                     kvas=[100.0, 20.0, 85.0] * um.kW,
    >>>                     taps=[1.0, 1.0, 1.0],
    >>>                     pctrs=[1.0, 1.0, 1.0]).aka('trf')
    >>> twc['a', 'b', 'c'] << trf
    >>> twc['b', 'bb'] << kp.Line(units= 'm', length=10.0 * um.m).aka('lbb')
    >>> twc['c', 'cc'] << kp.Line(units= 'm', length=15.0 * um.m).aka('lcc')
    >>> twc['bb',] << kp.Load(kv=2.0 * um.kV, kw = 5.0 * um.kW).aka('loadlow')
    >>> twc['cc',] << kp.Load(kv=3.0 * um.kV, kw = 35.0 * um.kW).aka('loadhi')
    >>> twc.solve()
    >>> bvo = kp.gv.VoltageView(twc)
    >>> print(round(bvo['bb'], 2))
    [ 974.96 -581.68j -991.23 -553.5j    16.27+1135.18j] volt
    >>> print(round(bvo['b','bb'], 2))
    [[ 977.18 -581.47j -992.16 -555.52j   14.98+1137.j  ] [ 974.96 -581.68j -991.23 -553.5j    16.27+1135.18j]] volt


.. _networkx: https://networkx.github.io/