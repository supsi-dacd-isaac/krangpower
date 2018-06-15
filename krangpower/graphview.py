import weakref
from functools import singledispatch

import networkx as nx

from .krangsuit import Krang

__all__ = ['GraphView']


class GraphView(nx.Graph):
    def __init__(self, busfun, edgefun, ckgr: Krang, raw_mode=False):
        """Applies a busfun and and edgefun to the graph of a Krang. It then is indicizable with bus names and tuples
        of two bus names to retrieve the relative results."""
        super().__init__(incoming_graph_data=None)

        wckgr = weakref.proxy(ckgr)

        gr = wckgr.graph()
        self.bus_pos = wckgr.bus_coords()
        stpos = {x: y for x, y in self.bus_pos.items() if y is not None}
        if stpos == {}:
            self.pad_pos = nx.spring_layout(gr)
        else:
            self.pad_pos = nx.spring_layout(gr, pos=stpos)

        self.raw_mode = raw_mode  # when false, getitem behavior is simplified
        """If True, __getitem__ has the same behavior as nx.Graph. If False, __getitem__ returns the values of busfun
        or edgefun directly."""

        if busfun is not None:
            self.bus_prop = busfun.__name__
            for n in gr.nodes:
                self.add_node(n, **{self.bus_prop: busfun(gr.nodes[n])})
        else:
            self.bus_prop = None
            for n in gr.nodes:
                self.add_node(n)

        if edgefun is not None:
            self.edge_prop = edgefun.__name__
            for e in gr.edges:
                self.add_edge(*e, **{self.edge_prop: edgefun(gr.edges[e])})
        else:
            self.edge_prop = None
            for e in gr.edges:
                self.add_edge(*e)

    def __getitem__(self, item):
        """Gets the value of busfun at a bus, or the value of edgefun at an edge.
        If GraphView.raw_mode == True, it behaves like a nx.Graph."""
        if self.raw_mode:
            return super().__getitem__(item)
        else:
            return _lean_getitem(item, self)


@singledispatch
def _lean_getitem(item, gr: GraphView):
    return gr.nodes[item][gr.bus_prop]


@_lean_getitem.register(tuple)
def _(item, gr: GraphView):
    return gr.edges[item][gr.edge_prop]
