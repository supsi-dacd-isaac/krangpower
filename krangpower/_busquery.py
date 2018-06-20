import numpy as np
from ._config_loader import UM


def get_fun(fun_name):
    g = globals()
    return g[fun_name]


def voltage(oek, busview, busname):
    return oek['bus.' + busname].Voltages()


def voltageangle(oek, busview, busname):
    return oek['bus.' + busname].VMagAngle()


def nloads(oek, busview, busname):
    return len([element for element in busview.content if element.type == 'load'])


def absvoltage(oek, busview, busname):
    return [np.abs(v) for v in voltage(oek, busview, busname)]


def totload(oek, busview, busname):
    els = [element for element in busview.content if element.type == 'load']
    tkw = 0.0 * UM.kW
    tkvar = 0.0 * UM.kVA
    for e in els:
        tkw += e['kw']
        tkvar += e['kvar']

    return tkw, tkvar

