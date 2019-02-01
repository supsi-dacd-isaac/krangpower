# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import weakref
from functools import singledispatch
from collections import OrderedDict

import networkx as nx

from ._krangsuit import Krang

__all__ = ['GraphView']


class GraphView(nx.Graph):
    def __init__(self, busfun, edgefun, ckgr: Krang, raw_mode=False):
        """Applies a busfun and and edgefun to the graph of a Krang. It then is indicizable with bus names and tuples
        of two bus names to retrieve the relative results."""
        super().__init__(incoming_graph_data=None)

        wckgr = weakref.proxy(ckgr)

        # self.gr = wckgr.graph()
        self._g_nodes = wckgr.graph().nodes
        self._g_edges = wckgr.graph().edges
        self.no = wckgr.brain.Circuit.AllBusNames()
        self.bus_pos = wckgr.bus_coords()
            
        self.raw_mode = raw_mode  # when false, getitem behavior is simplified
        """If True, __getitem__ has the same behavior as nx.Graph. If False, __getitem__ returns the values of busfun
        or edgefun directly."""

        if busfun is not None:
            self.bus_prop = busfun.__name__
            for n in self._g_nodes:
                self.add_node(n, **{self.bus_prop: busfun(self._g_nodes[n])})
        else:
            self.bus_prop = None
            for n in self._g_nodes:
                self.add_node(n)

        if edgefun is not None:
            self.edge_prop = edgefun.__name__
            for e in self._g_edges:
                self.add_edge(*e, **{self.edge_prop: edgefun(self._g_edges[e])})
        else:
            self.edge_prop = None
            for e in self._g_edges:
                self.add_edge(*e)

    def __getitem__(self, item):
        """Gets the value of busfun at a bus, or the value of edgefun at an edge.
        If GraphView.raw_mode == True, it behaves like a nx.Graph."""
        if self.raw_mode:
            return super().__getitem__(item)
        else:
            return _lean_getitem(item, self)
        
    # @property
    # def pad_pos(self):
    #     stpos = {x: y for x, y in self.bus_pos.items() if y is not None}
    #     if stpos == {}:
    #         return nx.spring_layout(self.gr)
    #     else:
    #         return nx.spring_layout(self.gr, pos=stpos)

    def get_edge_dict(self, convert_to_unit=None):
        if convert_to_unit is None:
            return {e: self[e] for e in self.edges}
        else:
            return {e: self[e].to(convert_to_unit).magnitude for e in self.edges}

    def get_node_dict(self, convert_to_unit=None):
        if convert_to_unit is None:
            ud = {n: self[n] for n in self.nodes}
            od = OrderedDict()
            for nodename in self.no:
                od[nodename] = ud[nodename]
            return od
        else:
            ud = {n: self[n].to(convert_to_unit).magnitude for n in self.nodes}
            od = OrderedDict()
            for nodename in self.no:
                od[nodename] = ud[nodename]
            return od


@singledispatch
def _lean_getitem(item, gr: GraphView):
    return gr.nodes[item][gr.bus_prop]


@_lean_getitem.register(tuple)
def _(item, gr: GraphView):
    return gr.edges[item][gr.edge_prop]
