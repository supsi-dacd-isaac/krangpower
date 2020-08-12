# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import weakref
from functools import singledispatch
from collections import OrderedDict

# import matplotlib
# matplotlib.use('Qt5Agg')
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.colors import Normalize
from matplotlib import pyplot as plt

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

    def plot(self,
             nodelist=None,
             edgelist=None,
             cmap_node='jet',
             cmap_edge='jet',
             **kwargs):

        stpos = {x: y for x, y in self.bus_pos.items() if y is not None}
        if stpos == {}:
            posi = nx.spring_layout(self)
        else:
            posi = nx.spring_layout(self, pos=stpos, fixed=stpos)

        if nodelist is None:
            nodelist = list(self.nodes)
        if edgelist is None:
            edgelist = list(self.edges)

        if self.bus_prop is not None:

            nd = self.get_node_dict()

            colors = [nd[x].magnitude for x in nodelist]
            unit = next(iter(nd.values())).units
            vmin = min(colors)
            vmax = max(colors)
            mcm = get_cmap(cmap_node)

            vmin = 242.45
            vmax= 242.55

            norm = Normalize(vmin=vmin, vmax=vmax)
            sm = ScalarMappable(norm=norm, cmap=cmap_node)
            sm.set_array(colors)

            nx.draw_networkx_nodes(self,
                                   pos=posi,
                                   nodelist=nodelist,
                                   node_color=colors,
                                   vmin=vmin, vmax=vmax,
                                   cmap=mcm,
                                   # style=(5, (5, 1)),
                                   **kwargs)

            cb = plt.colorbar(sm)
            cb.set_label(str(unit))
            # cb.set_array(colors)

        self._multiline_edge_plot(1, edgelist, posi)
        # plt.title(self.__class__.__name__)
        # plt.gca().axis('off')
        # plt.show(block=False)

    def _multiline_edge_plot(self, n, edgelist, posi, edge_color='k', base_color='w', base_width=0.8):
        fw = 2*n-1
        iswhite = [edge_color, base_color]
        for idl, ll in enumerate(range(fw, 0, -2)):
            nx.draw_networkx_edges(self,
                                   pos=posi,
                                   edgelist=edgelist,
                                   width=ll*base_width,
                                   edge_color=iswhite[idl % 2])


@singledispatch
def _lean_getitem(item, gr: GraphView):
    return gr.nodes[item][gr.bus_prop]


@_lean_getitem.register(tuple)
def _(item, gr: GraphView):
    return gr.edges[item][gr.edge_prop]
