# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import numpy as np

from .._graphview import GraphView
from .._krangsuit import Krang


__all__ = ['BusVoltageView', 'VoltageView', 'CurrentView', 'BaseVoltageView']


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

    def min_v(self):
        vs = [self[x] for x in self.nodes]
        return min(vs)


class CurrentView(GraphView):
    def __init__(self, ckgr: Krang):

        def edgeI(edg):
            return edg['el'][0].Currents()

        super().__init__(None, edgeI, ckgr)


class BaseVoltageView(GraphView):
    def __init__(self, ckgr: Krang, voltage_bases):

        ckgr.set(voltagebases=voltage_bases)
        ckgr.command('calcvoltagebases')

        def basebuskV(bus):
            return bus['bus'].kVBase()

        super().__init__(basebuskV, None, ckgr)
