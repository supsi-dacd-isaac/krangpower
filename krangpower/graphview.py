import networkx as nx
from functools import singledispatch
import numpy as np
from krangpower import Krang

__all__ = ['GraphView', 'CurrentView', 'VoltageView', 'BusVoltageView']


class GraphView(nx.Graph):
    def __init__(self, busfun, edgefun, ckgr: Krang):
        super().__init__(incoming_graph_data=None)

        gr = ckgr.graph()
        self.bus_pos = ckgr.bus_coords()
        stpos = {x: y for x, y in self.bus_pos.items() if y is not None}
        if stpos == {}:
            self.pad_pos = nx.spring_layout(gr)
        else:
            self.pad_pos = nx.spring_layout(gr, pos=stpos)

        self.raw_mode = False  # when false, getitem behavior is simplified

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


class VoltageView(GraphView):
    def __init__(self, ckgr: Krang):
        def busV(bus):
            return bus['bus'].Voltages()

        def edgeV(edg):
            return edg['el'][0].Voltages()

        super().__init__(busV, edgeV, ckgr)


class BusVoltageView(GraphView):
    def __init__(self, ckgr: Krang):
        def avgbusV(bus):
            return np.mean(np.abs(bus['bus'].Voltages()))

        super().__init__(avgbusV, None, ckgr)

    def xyz(self):
        return np.asarray([x[0] for x in self.pad_pos.values()]), \
               np.asarray([x[1] for x in self.pad_pos.values()]), \
               np.asarray([self[k].magnitude for k in self.pad_pos.keys()])


class CurrentView(GraphView):
    def __init__(self, ckgr: Krang):

        def edgeI(edg):
            return edg['el'][0].Currents()

        super().__init__(None, edgeI, ckgr)
