import numpy as np

from krangpower.graphview import GraphView
from krangpower.krangsuit import Krang


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