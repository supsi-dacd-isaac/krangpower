# OpendssdirectEnhancer by Federico Rosato
# a wrapper for opendssdirect.py by Dheepak Krishnamurthy and Maximilian J. Zangs

from functools import reduce
import copy
from math import sqrt
from operator import getitem
from typing import Callable

import numpy as np
import opendssdirect as odr
from pandas import DataFrame
import components as co
from components import um, logger
from components import _resolve_unit, _SnpMatrix, _pint_qty_type, _odssrep, _type_recovery

from utils.aux_fcn import lower as _lower
from utils.aux_fcn import pairwise as _pairwise

__all__ = ['OpendssdirectEnhancer']

# the default entity parameter values are loaded in order to allow correct type casting and comparison with the default
# when calling _Packed's __getitem__ method
_default_entities = co.default_comp


# <editor-fold desc="Auxiliary functions">

def _assign_unit(item, unit: type(um.m) or None):
    if unit is None:
        return item
    elif isinstance(item, dict):
        return {k: v * unit for k, v in item.items()}
    elif isinstance(item, DataFrame):
        # pandas' dataframe is a mess together with pint
        return item
    # elif hasattr(item, '__iter__'):
    #     # return [el * unit for el in item]
    #     return _asarray(item) * unit
    else:
        return item * unit


def _asarray(item):
    return np.asarray(item)


def _couplearray(item):
    return np.array(item[0::2]) + np.array(item[1::2]) * 1j


def _terminalize(item):

    # note that, when I pass an item to terminalize, I am taking for granted that I can find nterm and ncond in the
    # respective calls to odr. If you called odr.CktElement.Powers(), for example, I take it that you knew what
    # you were doing. Calls coming from PackedElements, instead, should perform the cktelement selection just before
    # the call.

    nterm = odr.CktElement.NumTerminals()
    ncond = odr.CktElement.NumConductors()

    assert len(item) == nterm * ncond * 2
    cpxr = np.zeros([nterm, ncond], 'complex')

    for idx, couple in enumerate(_pairwise(item)):
        real = couple[0]
        imag = couple[1]
        cpxr[int(idx / ncond), (idx % ncond)] = np.sum([np.multiply(1j, imag), real], axis=0)
        cpxr[int(idx / ncond), (idx % ncond)] = np.sum([np.multiply(1j, imag), real], axis=0)

    return cpxr


def _cpx(item):
    return item[0] + 1j*item[1]


def _dictionize_cktbus(item):
    return dict(zip(_lower(odr.Circuit.AllBusNames()), item))


def _dictionize_cktels(item):
    return dict(zip(_lower(odr.Circuit.AllElementNames()), item))


def _dictionize_cktnodes(item):
    return dict(zip(_lower(odr.Circuit.YNodeOrder()), item))


def _matricize_ybus(item):

    raw_n_ord = _lower(odr.Circuit.YNodeOrder())

    mtx = _matricize(item)

    return DataFrame(data=mtx, index=raw_n_ord, columns=raw_n_ord)


def _matricize(item):

    if len(item) == 1:
        return item

    side = sqrt(len(item)/2)
    assert side == int(side)
    side = int(side)
    mtx = np.reshape(item[0::2], (side, side)) + \
        np.reshape(item[1::2], (side, side)) * 1j

    return np.transpose(mtx)  # it's transposed because it's originally given in column order


def _cast_dumbstring(string: str, type):

    if type == str:
        return string
    if type in (int, float):
        return type(string)
    elif type == np.matrix:
        return _SnpMatrix(string
                          .replace(' |', ';')
                          .replace('|', ';')
                          .replace('[', '')
                          .replace(' ]', '')
                          .replace(']', '')
                          .replace(' ', ','))
    else:
        raise TypeError('Could not cast the DSS property string: type {0} unknown'.format(str(type)))


# there's no XYCurves.AllNames() or similar, so we have to mock up one ourselves
def _xycurve_names():
    xynames = []
    i = 1
    odr.XYCurves.First()
    while i != 0:
        xynames.append(odr.XYCurves.Name())
        i = odr.XYCurves.Next()

    return xynames
# </editor-fold>


class OpendssdirectEnhancer:
    """This class is designed to behave exactly as the opendssdirect MODULE, but adds some handy functionality:

    -   Items returned as a list of floats (e.g., <>.Circuit.Losses(), <>.Bus.Voltages()...)
        are returned as lists of complex numbers, matrices of [nterm x ncond], etc. as appropriate
    -   Structured items such as opendssdirect.Circuit.SistemY() are returned as pandas.DataFrame
    -   Items come, where appropriate, as quantities (from the pint package) with the appropriate measurement unit

    -   Supports getitem with the fully qualified names of the circuit's elements in the way shown below:

            >>> MOE = OpendssdirectEnhancer()
            >>> MOE.utils.run_command("Clear")
            ''
            >>> MOE.utils.run_command("New object = circuit.myckt bus1=sourcebus basekv=11.0 pu=1.0 angle=0.0 phases=3")
            ''
            >>> MOE.utils.run_command("New load.myload bus1=sourcebus kw=10.0 kv=11.0 basefreq=50.0")
            ''
            >>> MOE.utils.run_command("solve")
            ''
            >>> MOE['load.myload'].CFactor() # taken from the Loads interface
            4.0
            >>> abs(MOE['load.myload'].Currents()[0,1]) # taken from the CktElement interface
            <Quantity(0.5964385133615372, 'ampere')>

    Note that the user needs not to worry about selecting the item with MOE.Circuit.SetActiveElement('load.myload'),
    MOE.Laods.Name('load.myload'); the enhancer takes care of the selections automatically.
    """

    line_um = {
        0: um.unitlength,
        1: um.mile,
        2: um.kft,
        3: um.km,
        4: um.m,
        5: um.ft,
        6: um.inch,
        7: um.cm
    }

    trt = {'Settings':
                {'VoltageBases': (_asarray,)},
           'Bus':
                {'CplxSeqVoltages': (_couplearray,),
                 'SeqVoltages': (_couplearray,),
                 'Voltages': (_couplearray,),
                 'YscMatrix': (_asarray, _matricize,),
                 'Zsc0': (_cpx,),
                 'Zsc1': (_cpx,),
                 'ZscMatrix': (_asarray, _matricize,),
                 # puvmagangle
                 },

           'Circuit':
                {'AllBusDistances': (_dictionize_cktbus,),
                 'AllBusVMag': (_dictionize_cktbus,),
                 'AllBusVolts': (_couplearray, _dictionize_cktbus),
                 'AllElementLosses': (_couplearray, _dictionize_cktels),
                 'LineLosses': (_cpx,),
                 'Losses': (_cpx,),
                 'TotalPower': (_cpx,),
                 'YCurrents': (_couplearray, _dictionize_cktnodes),
                 'YNodeVArray': (_couplearray, _dictionize_cktnodes),
                 'SystemY': (_matricize_ybus,)},

           'CktElement':
                {'CplxSeqCurrents': (_couplearray,),
                 'CplxSeqVoltages': (_couplearray,),
                 'Currents': (_terminalize,),
                 'CurrentsMagAng': (_terminalize,),
                 'Losses': (_cpx,),
                 'PhaseLosses': (_couplearray,),
                 'Powers': (_terminalize,),
                 'SeqCurrents': (_couplearray,),
                 'SeqVoltages': (_couplearray,),
                 'Voltages': (_terminalize,),
                 'VoltagesMagAng': (_terminalize,),
                 'Yprim': (_matricize,)
                 }
        # todo complete with the other items
    }

    umr = {'Settings':
                {'VoltageBases': um.kV},
           'Bus':
                {'CplxSeqVoltages': um.V,
                 'Cust_Duration': um.hr,
                 'Distance': um.km,
                 'Int_Duration': um.hr,
                 'Isc': um.A,
                 'SeqVoltages': um.V,
                 'TotalMiles': um.mile,
                 'VLL': um.V,
                 'VMagAngle': um.deg,
                 'Voc': um.V,
                 'Voltages': um.V,
                 'YscMatrix': um.siemens,
                 'Zsc0': um.ohm,
                 'Zsc1': um.ohm,
                 'ZscMatrix': um.ohm,
                 'kVBase': um.kV
                 #puvmagangle
                 },
           'CapControls':
                {'Vmin': um.V,
                 'Vmax': um.V},
           'Capacitors':
                {'kv': um.kV,
                 'kvar': um.kVA},
           'Circuit':
                {'AllBusDistances': um.km,
                 'AllBusVMag': um.V,
                 'AllBusVolts': um.V,
                 'AllElementLosses': um.kW,
                 'LineLosses': um.W,
                 'Losses': um.W,
                 'TotalPower': um.kVA,
                 'YCurrents': um.A,
                 'YNodeVArray': um.V,
                 'SystemY': um.siemens},
           'Fuses':
                {'RatedCurrent': um.A},
           'Generators':
                {'kV': um.kV,
                 'kVARated': um.kVA,
                 'kW': um.kW,
                 'kvar': um.kVA},
           'Isource':
                {'Amps': um.A,
                 'AngleDeg': um.deg,
                 'Frequency': um.Hz},
           'Loads':
                {'PctMean': um.pct,
                 'PctStdDev': um.pct,
                 'Rneut': um.ohm,
                 'XfkVA': um.kVA,
                 'Xneut': um.ohm,
                 'kV': um.kVA,
                 'kVABase': um.kVA,
                 'kW': um.kW,
                 'kWh': um.kWh,
                 'kWhDays': um.day,
                 'kvar': um.kVA,
                 'puSeriesRL': um.pct},
           'Meters':
                {
                 'AvgRepairTime': um.hr,
                 'CalcCurrent': um.A,
                 'PeakCurrent': um.A},
           # todo monitor bytestream...
           'PDElements':
                {'PctPermanent': um.pct,
                 'RepairTime': um.hr,
                 'TotalMiles': um.mile},
           'PVSystems':
                {'Irradiance': um.W / (um.m ** 2),
                 'kVARated': um.kVA,
                 'kW': um.kW,
                 'kvar': um.kVA},
           'Solution':
                {'Capkvar': um.kVA,
                 'DblHour': um.hr,
                 'Frequency': um.Hz,
                 'GenkW': um.kW,
                 'Hour': um.hr,
                 'Seconds': um.s,
                 'StepSize': um.s,
                 'StepSizeHr': um.hr,
                 'StepSizeMin': um.min,
                 'TimeTimeStep': um.microsecond,
                 'TotalTime': um.microsecond},
           'CktElement':
                {'CplxSeqCurrents': um.A,
                 'CplxSeqVoltages': um.V,
                 'Currents': um.A,
                 'CurrentsMagAng': um.deg,
                 'Emergamps': um.A,
                 'Losses': um.W,
                 'NormalAmps': um.A,
                 'PhaseLosses': um.kW,
                 'Powers': um.kW,
                 'SeqCurrents': um.A,
                 'SeqVoltages': um.V,
                 'Voltages': um.V,
                 'VoltagesMagAng': um.deg,
                 'Yprim': um.siemens
                 }
        # todo finish the other items
    }

    def __init__(self, stack=None):
        if stack is None:
            self.stack = []
        else:
            self.stack = stack

    def __getattr__(self, item):
        if self._chain_hasattr(self.stack + [item]):
            return OpendssdirectEnhancer(self.stack + [item])
            # return getattr(odr, item)
        else:
            raise AttributeError('Could not find the attribute {0} in the available interfaces.'.format(item))

    def __getitem__(self, item):
        """bracket indicization looks for an object with the name desired and returns a nice packedelement that you
        can call with all the usual attributes, but without ever worrying again about SetActiveElement() and the likes.
        """
        if item.lower() not in map(lambda name: name.lower(), self.get_all_names()):
            raise KeyError('Element {0} was not found in the circuit'.format(item))

        return _PackedOpendssElement(*item.lower().split('.', 1), self)

    def __str__(self):
        return '<OpendssdirectEnhancer -> {0}>'.format('.'.join(['opendssdirect'] + self.stack))

    def __repr__(self):
        return self.__str__()

    def __call__(self, *args):

        if self.stack[0] == 'Lines':
            um_d = self.line_umd(self.line_um[odr.Lines.Units()])
        elif self.stack[0] == 'LoadShape':
            um_d = self.loadshape_umd(self.line_um[bool(odr.LoadShape.UseActual())])
        else:
            um_d = self.umr

        try:
            ums = reduce(getitem, self.stack, um_d)  # nested dict search
        except KeyError:
            ums = None

        try:
            trt = reduce(getitem, self.stack, self.trt)
        except KeyError:
            trt = tuple()  # empty tuple in order to shortcut iteration

        odrobj = self._chain_getattr()
        e_ordobj = odrobj(*args)
        for t in trt:
            e_ordobj = t(e_ordobj)

        return _assign_unit(e_ordobj, ums)

    def _jsonize(self):
        pass

    @staticmethod
    def get_all_names():

        anl = []
        anl.extend(map(lambda bn: 'bus.' + bn, odr.Circuit.AllBusNames()))
        anl.extend(odr.Circuit.AllElementNames())
        anl.extend(map(lambda ln: 'loadshape.' + ln, odr.LoadShape.AllNames()))
        anl.extend(map(lambda ln: 'xycurve.' + ln, _xycurve_names()))

        return anl

    @staticmethod
    def loadshape_umd(use_actual: bool):

        if use_actual:
            return {'Loadshapes':
                        {'HrInterval': um.hr,
                         'MinInterval': um.min,
                         'PBase': um.kW,
                         'PMult': um.kW,
                         'QBase': um.kVA,
                         'QMult': um.kVA,
                         'SInterval': um.s,
                         'TimeArray': um.hr}}

        else:
            return {'Loadshapes':
                        {'HrInterval': um.hr,
                         'MinInterval': um.min,
                         'PBase': um.kW,
                         'PMult': um.dimensionless,
                         'QBase': um.kW,
                         'QMult': um.dimensionless,
                         'SInterval': um.s,
                         'TimeArray': um.hr}}

    @staticmethod
    def line_umd(unit_length):
        """Special case of call if the object is a line; needed because it's a type of object for which a "units"
        properties exists that influences the other quantities dimensions."""

        line_qties = {'Lines': {'C0': um.nF/unit_length,
                                'C1': um.nF/unit_length,
                                'CMatrix': um.nF/unit_length,
                                'EmergAmps': um.A,
                                'Length': unit_length,
                                'NormAmps': um.A,
                                'R0': um.ohm/unit_length,
                                'R1': um.ohm/unit_length,
                                'RMatrix': um.ohm/unit_length,
                                'Rg': um.ohm/unit_length,
                                'Rho': um.ohm/unit_length,
                                'X0': um.ohm/unit_length,
                                'X1': um.ohm/unit_length,
                                'XMatrix': um.ohm/unit_length,
                                'Xg': um.ohm/unit_length,
                                'Yprim': um.siemens/unit_length}}

        return line_qties

    def _chain_getattr(self, stack=None, start_obj=odr):

        if stack is None:
            stack = self.stack

        if len(stack) == 1:
            return getattr(start_obj, stack[0])
        else:
            return self._chain_getattr(stack[1:], getattr(start_obj, stack[0]))

    def _chain_hasattr(self, stack=None, start_obj=odr):

        if stack is None:
            stack = self.stack

        if len(stack) == 1:
            return hasattr(start_obj, stack[0])
        else:
            return self._chain_getattr(stack[1:], getattr(start_obj, stack[0]))


# p_odr = odr

class _PackedOpendssElement:

    def __init__(self, eltype, name, p_odr):

        # todo complete and verify these dicts
        interfaces = {'bus': (p_odr.Bus,),
                      'load': (p_odr.Loads, p_odr.CktElement),
                      'monitor': (p_odr.Monitors, p_odr.CktElement),
                      'line': (p_odr.Lines, p_odr.CktElement, p_odr.PDElements),
                      'capcontrol': (p_odr.CapControls,),
                      'capacitor': (p_odr.Capacitors, p_odr.CktElement, p_odr.PDElements),
                      'isource': (p_odr.Isource, p_odr.CktElement),
                      'meter': (p_odr.Meters, p_odr.CktElement),
                      'vsource': (p_odr.Vsources, p_odr.CktElement),
                      'pvsystem': (p_odr.PVsystems, p_odr.CktElement),
                      'regcontrol': (p_odr.RegControls,),
                      'xycurve': (p_odr.XYCurves,),
                      'transformer': (p_odr.Transformers, p_odr.CktElement, p_odr.PDElements),
                      'storage': (p_odr.CktElement,),
                      'loadshape': (p_odr.LoadShape,),
                      'reactor': (p_odr.PDElements,)
                      }

        selectors = {'load': (p_odr.Loads.Name, p_odr.Circuit.SetActiveElement),
                     'bus': (p_odr.Circuit.SetActiveBus,),
                     'monitor': (p_odr.Monitors.Name, p_odr.Circuit.SetActiveElement),
                     'line': (p_odr.Lines.Name, p_odr.Circuit.SetActiveElement, p_odr.PDElements.Name),
                     'capcontrol': (p_odr.CapControls,),
                     'capacitor': (p_odr.Capacitors.Name, p_odr.Circuit.SetActiveElement, p_odr.PDElements.Name),
                     'isource': (p_odr.Isource.Name, p_odr.Circuit.SetActiveElement),
                     'meter': (p_odr.Meters.Name,),
                     'vsource': (p_odr.Vsources.Name, p_odr.Circuit.SetActiveElement),
                     'pvsystem': (p_odr.PVsystems.Name, p_odr.Circuit.SetActiveElement),
                     'regcontrol': (p_odr.RegControls.Name,),
                     'xycurve': (p_odr.XYCurves.Name,),
                     'transformer': (p_odr.Transformers.Name(), p_odr.Circuit.SetActiveElement),
                     'storage': (p_odr.Circuit.SetActiveElement,),
                     'loadshape': (p_odr.LoadShape.Name()),
                     'reactor': (p_odr.PDElements.Name,)}

        interface_methods = {('ActiveClass',): ['ActiveClassName', 'AllNames', 'Count', 'First', 'Name', 'Next', 'NumElements'],
                             ('Basic',): ['AllowForms', 'Classes', 'ClearAll', 'DataPath', 'DefaultEditor', 'NewCircuit', 'NumCircuits', 'NumClasses', 'NumUserClasses', 'Reset', 'ShowPanel', 'Start', 'UserClasses', 'Version'],
                             ('Bus',): ['Coorddefined', 'CplxSeqVoltages', 'Distance', 'GetUniqueNodeNumber', 'Isc', 'Lambda', 'Name', 'Nodes', 'NumNodes', 'PuVoltage', 'SectionID', 'SeqVoltages', 'TotalMiles', 'VLL', 'VMagAngle', 'Voc', 'Voltages', 'X', 'Y', 'YscMatrix', 'Zsc0', 'Zsc1', 'ZscMatrix', 'ZscRefresh', 'kVBase', 'puVLL', 'puVmagAngle'],
                             ('Capacitors',): ['AddStep', 'AllNames', 'AvailableSteps', 'Close', 'Count', 'First', 'IsDelta', 'Name', 'Next', 'NumSteps', 'Open', 'States', 'SubtractStep', 'kV', 'kvar'],
                             ('CapControls',): ['AllNames', 'CTRatio', 'Capacitor', 'Count', 'Delay', 'DelayOff', 'First', 'Mode', 'MonitoredObj', 'MonitoredTerm', 'Name', 'Next', 'OFFSetting', 'ONSetting', 'PTRatio', 'UseVoltOverride', 'Vmax', 'Vmin'],
                             ('Circuit',): ['AllBusDistances', 'AllBusMagPu', 'AllBusNames', 'AllBusVMag', 'AllBusVolts', 'AllElementLosses', 'AllElementNames', 'AllNodeDistances', 'AllNodeNames', 'Capacity', 'Disable', 'Enable', 'EndOfTimeStepUpdate', 'FirstElement', 'FirstPCElement', 'FirstPDElement', 'LineLosses', 'Losses', 'Name', 'NextElement', 'NextPCElement', 'NextPDElement', 'NumBuses', 'NumCktElements', 'NumNodes', 'ParentPDElement', 'Sample', 'SaveSample', 'SetActiveBus', 'SetActiveBusi', 'SetActiveClass', 'SetActiveElement', 'SubstationLosses', 'TotalPower', 'UpdateStorage', 'YCurrents', 'YNodeOrder', 'YNodeVArray'],
                             ('CktElement',): ['AllPropertyNames', 'AllVariableNames', 'AllVariableValues', 'BusNames', 'Close', 'CplxSeqCurrents', 'CplxSeqVoltages', 'Currents', 'CurrentsMagAng', 'DisplayName', 'EmergAmps', 'Enabled', 'EnergyMeter', 'GUID', 'HasSwitchControl', 'HasVoltControl', 'IsOpen', 'Losses', 'Name', 'NodeOrder', 'NormalAmps', 'NumConductors', 'NumControls', 'NumPhases', 'NumProperties', 'NumTerminals', 'OCPDevIndex', 'OCPDevType', 'Open', 'PhaseLosses', 'Powers', 'Residuals', 'SeqCurrents', 'SeqPowers', 'SeqVoltages', 'Variablei', 'Voltages', 'VoltagesMagAng', 'YPrim'],
                             ('dss',): ['ActiveClass', 'Basic', 'Bus', 'CapControls', 'Capacitors', 'Circuit', 'CktElement', 'Element', 'Executive', 'Fuses', 'Generators', 'Isource', 'Lines', 'LoadShape', 'Loads', 'Meters', 'Monitors', 'PDElements', 'PVsystems', 'Parser', 'Properties', 'Reclosers', 'RegControls', 'Relays', 'Sensors', 'Settings', 'Solution', 'SwtControls', 'Topology', 'Transformers', 'Vsources', 'XYCurves', 'dss'],
                             ('Element',): ['AllPropertyNames', 'Name', 'NumProperties'],
                             ('Executive',): ['Command', 'CommandHelp', 'NumCommands', 'NumOptions', 'Option', 'OptionHelp', 'OptionValue'],
                             ('Fuses',): ['AllNames', 'Close', 'Count', 'First', 'Idx', 'IsBlown', 'MonitoredObj', 'MonitoredTerm', 'Name', 'Next', 'NumPhases', 'Open', 'RatedCurrent', 'SwitchedObj', 'TCCCurve'],
                             ('Generators',): ['AllNames', 'Count', 'First', 'ForcedON', 'Idx', 'Model', 'Name', 'Next', 'PF', 'Phases', 'RegisterNames', 'RegisterValues', 'Vmaxpu', 'Vminpu', 'kV', 'kVARated', 'kW', 'kvar'],
                             ('Isource',): ['AllNames', 'Amps', 'AngleDeg', 'Count', 'First', 'Frequency', 'Name', 'Next'],
                             ('Lines',): ['AllNames', 'Bus1', 'Bus2', 'C0', 'C1', 'CMatrix', 'Count', 'EmergAmps', 'First', 'Geometry', 'Length', 'LineCode', 'Name', 'Next', 'NormAmps', 'NumCust', 'Parent', 'Phases', 'R0', 'R1', 'RMatrix', 'Rg', 'Rho', 'Spacing', 'Units', 'X0', 'X1', 'XMatrix', 'Xg', 'Yprim'],
                             ('Loads',): ['AllNames', 'AllocationFactor', 'CFactor', 'CVRCurve', 'CVRvars', 'CVRwatts', 'Class', 'Count', 'Daily', 'Duty', 'First', 'Growth', 'Idx', 'IsDelta', 'Model', 'Name', 'Next', 'NumCust', 'PF', 'PctMean', 'PctStdDev', 'RelWeighting', 'Rneut', 'Spectrum', 'Status', 'Vmaxpu', 'VminEmerg', 'VminNorm', 'Vminpu', 'XfkVA', 'Xneut', 'Yearly', 'ZipV', 'kV', 'kVABase', 'kW', 'kWh', 'kWhDays', 'kvar', 'puSeriesRL'],
                             ('LoadShape',): ['AllNames', 'Count', 'First', 'HrInterval', 'MinInterval', 'Name', 'Next', 'Normalize', 'Npts', 'PBase', 'PMult', 'QBase', 'QMult', 'SInterval', 'TimeArray', 'UseActual'],
                             ('Meters',): ['AllBranchesInZone', 'AllEndElements', 'AllNames', 'AllocFactors', 'AvgRepairTime', 'CalcCurrent', 'CloseAllDIFiles', 'Count', 'CountBranches', 'CountEndElements', 'CustInterrupts', 'DIFilesAreOpen', 'DoReliabilityCalc', 'FaultRateXRepairHrs', 'First', 'MeteredElement', 'MeteredTerminal', 'Name', 'Next', 'NumSectionBranches', 'NumSectionCustomers', 'NumSections', 'OCPDeviceType', 'OpenAllDIFiles', 'PeakCurrent', 'RegisterNames', 'RegisterValues', 'Reset', 'ResetAll', 'SAIDI', 'SAIFI', 'SAIFIkW', 'Sample', 'SampleAll', 'Save', 'SaveAll', 'SectSeqidx', 'SectTotalCust', 'SeqListSize', 'SequenceList', 'SetActiveSection', 'SumBranchFltRates', 'TotalCustomers', 'Totals'],
                             ('Monitors',): ['AllNames', 'ByteStream', 'Count', 'Element', 'FileName', 'FileVersion', 'First', 'Mode', 'Name', 'Next', 'Process', 'ProcessAll', 'Reset', 'ResetAll', 'Sample', 'SampleAll', 'Save', 'SaveAll', 'Show', 'Terminal'],
                             ('PDElements',): ['AccumulatedL', 'Count', 'FaultRate', 'First', 'FromTerminal', 'IsShunt', 'Lambda', 'Name', 'Next', 'NumCustomers', 'ParentPDElement', 'PctPermanent', 'RepairTime', 'SectionID', 'TotalCustomers', 'TotalMiles'],
                             ('Properties',): ['Description', 'Name', 'Value'],
                             ('PVsystems',): ['AllNames', 'Count', 'First', 'Idx', 'Irradiance', 'Name', 'Next', 'kVARated', 'kW', 'kvar', 'pf'],
                             ('Reclosers',): ['AllNames', 'Close', 'Count', 'First', 'GroundInst', 'GroundTrip', 'Idx', 'MonitoredObj', 'MonitoredTerm', 'Name', 'Next', 'NumFast', 'Open', 'PhaseInst', 'PhaseTrip', 'RecloseIntervals', 'Shots', 'SwitchedObj', 'SwitchedTerm'],
                             ('RegControls',): ['AllNames', 'CTPrimary', 'Count', 'Delay', 'First', 'ForwardBand', 'ForwardR', 'ForwardVreg', 'ForwardX', 'IsInverseTime', 'IsReversible', 'MaxTapChange', 'MonitoredBus', 'Name', 'Next', 'PTRatio', 'ReverseBand', 'ReverseR', 'ReverseVreg', 'ReverseX', 'TapDelay', 'TapNumber', 'TapWinding', 'Transformer', 'VoltageLimit', 'Winding'],
                             ('Relays',): ['AllNames', 'Count', 'First', 'Idx', 'MonitoredObj', 'MonitoredTerm', 'Name', 'Next', 'SwitchedObj', 'SwitchedTerm'],
                             ('Sensors',): ['AllNames', 'Count', 'Currents', 'First', 'IsDelta', 'MeteredElement', 'MeteredTerminal', 'Name', 'Next', 'PctError', 'Reset', 'ResetAll', 'ReverseDelta', 'Weight', 'kVBase', 'kW', 'kvar'],
                             ('Settings',): ['AllocationFactors', 'AllowDuplicates', 'AutoBusList', 'CktModel', 'EmergVmaxpu', 'EmergVminpu', 'LossRegs', 'LossWeight', 'NormVmaxpu', 'NormVminpu', 'PriceCurve', 'PriceSignal', 'Trapezoidal', 'UERegs', 'UEWeight', 'VoltageBases', 'ZoneLock'],
                             ('Solution',): ['AddType', 'Algorithm', 'BuildYMatrix', 'Capkvar', 'CheckControls', 'CheckFaultStatus', 'Cleanup', 'ControlActionsDone', 'ControlIterations', 'ControlMode', 'Converged', 'Convergence', 'DblHour', 'DefaultDaily', 'DefaultYearly', 'DoControlActions', 'EventLog', 'FinishTimeStep', 'Frequency', 'GenMult', 'GenPF', 'GenkW', 'Hour', 'InitSnap', 'Iterations', 'LDCurve', 'LoadModel', 'LoadMult', 'MaxControlIterations', 'MaxIterations', 'Mode', 'ModeID', 'MostIterationsDone', 'Number', 'PctGrowth', 'ProcessTime', 'Random', 'SampleControlDevices', 'SampleDoControlActions', 'Seconds', 'Solve', 'SolveDirect', 'SolveNoControl', 'SolvePFlow', 'SolvePlusControl', 'StepSize', 'StepSizeHr', 'StepSizeMin', 'SystemYChanged', 'TimeTimeStep', 'TotalIterations', 'TotalTime', 'Year'],
                             ('SwtControls',): ['Action', 'AllNames', 'Count', 'Delay', 'First', 'IsLocked', 'Name', 'Next', 'SwitchedObj', 'SwitchedTerm'],
                             ('Topology',): ['ActiveBranch', 'ActiveLevel', 'AllIsolatedBranches', 'AllIsolatedLoads', 'AllLoopedPairs', 'BranchName', 'BusName', 'First', 'FirstLoad', 'ForwardBranch', 'LoopedBranch', 'Next', 'NextLoad', 'NumIsolatedBranches', 'NumIsolatedLoads', 'NumLoops', 'ParallelBranch'],
                             ('Transformers',): ['AllNames', 'Count', 'First', 'IsDelta', 'MaxTap', 'MinTap', 'Name', 'Next', 'NumTaps', 'NumWindings', 'R', 'Rneut', 'Tap', 'Wdg', 'XfmrCode', 'Xhl', 'Xht', 'Xlt', 'Xneut', 'kV', 'kVA'],
                             ('utils',): ['Iterator'],
                             ('Vsources',): ['AllNames', 'AngleDeg', 'BasekV', 'Count', 'First', 'Frequency', 'Name', 'Next', 'PU', 'Phases'],
                             ('XYCurves',): ['Count', 'First', 'Name', 'Next', 'Npts', 'X', 'XArray', 'XScale', 'XShift', 'Y', 'YArray', 'YScale', 'YShift']}

        self._available_interfaces = interfaces[eltype]
        self._name = name
        # full qualified name
        # todo perform sanity checks on name
        self._selectors = selectors[eltype]
        # todo distinct linecode_a, etc and make sure that _eltype is coherent with components
        self._eltype = eltype
        self._p_odr = p_odr

        # dynamically add methods to expose
        for i in self._available_interfaces:
            for m in interface_methods.get(tuple(i.__name__.split('.')[-1]), []):
                setattr(self, m, self._craft_member(m))

    @property
    def topological(self):
        top_par_names = _default_entities['default_' + self._eltype]['topological']
        return tuple(self[t] for t in top_par_names.keys())

    @property
    def type(self):
        return self._eltype

    @property
    def fullname(self):
        return self._eltype + '.' + self._name

    def _craft_member(self, item: str):
        
        def pckdss_mthd(*args):
            # self._side_getattr(m) is a CallFinalizer
            return self._side_getattr(item)(*args)

        # poo_mthd.__repr__ = lambda: '<function {0} of {1}.{2}>'.format(item, self._eltype, self._name)
        # poo_mthd.__name__ = lambda: item
        return pckdss_mthd

    def _side_getattr(self, item):

        for itf in self._available_interfaces:
            if hasattr(itf, item):
                return _CallFinalizer(getattr(itf, item), self._selectors, self._name, str(itf))
                # break unnecessary
            else:
                continue
        else:
            raise AttributeError('"{0}" is neither an attribute of the _PackedOpendssElement nor of the {1}-type object'
                                 ' it wraps.'
                                 .format(item, self._eltype.upper()))

    def dump(self):
        """Returns a dict with all the properties-values pairs of the object as they would be returned
        by __getitem__."""
        try:
            props = self._p_odr[self.fullname].AllPropertyNames()
        except AttributeError:
            raise ValueError('{0}-type objects are not dumpable.'.format(self._eltype.upper()))
        return {p: self._p_odr[self.fullname][p] for p in props}

    def unpack(self, verbose=False):

        # todo use associated

        classmap = co.get_classmap()
        myclass = classmap[self._eltype]

        all_props = self.dump()
        valid_props = copy.deepcopy(_default_entities['default_' + self._eltype]['properties'])

        valid_props.update(_default_entities['default_' + self._eltype].get('associated', {}))

        if verbose:
            dep_prop = {k.lower(): v for k, v in all_props.items() if k.lower() in valid_props.keys()}
        else:
            dep_prop = {k.lower(): v for k, v in all_props.items() if k.lower() in valid_props.keys() and np.matrix(v != valid_props[k.lower()]).any()}

        if myclass.isnamed():
            obj = myclass(self._name, **dep_prop)
        else:
            obj = myclass(**dep_prop)

        obj.name = self._name

        return obj

    def __getattr__(self, item):
        return self._side_getattr(item)

    def __getitem__(self, item):
        """Gets a property of the item from the text interface (as opposed to a natively exposed attribute of the
        object's interface). For example, <>['xrharm'], where <> is a _PackedOpendssElement representing a Load,
        lets you access the 'xrharm' property, which is not available in any way through the standard interface."""

        rslt = self._p_odr.utils.run_command('? {0}.{1}'.format(self.fullname, item))

        if rslt == 'Property Unknown':  # this is what odr returns, instead of raising, if you request invalid props
            raise KeyError('Invalid property "{0}" for the object "{1}"'.format(item, self.fullname))

        # ok, so if we're here, rslt is a true value. Rslt, being retrieved by the text interface,
        # is now a raw string. We have to realize from the default entities its type and possible unit
        # and then act accordingly

        try:
            data_type = self._get_datatype(item)
        except KeyError:
            # happens, for example, when topological parameters - excluded from the default dicts - are requested
            return rslt

        rslt = _cast_dumbstring(rslt, data_type)

        unt = self._get_builtin_units(item)
        if unt is None:
            return rslt
        else:
            return rslt * unt

    def _get_datatype(self, item):
        return type(_default_entities['default_' + self._eltype]['properties'][item.lower()])

    def _get_builtin_units(self, item):
        raw_unit = _default_entities['default_' + self._eltype]['units'].get(item.lower(), None)
        if raw_unit is None:
            return None
        else:
            return _resolve_unit(raw_unit, self._get_matching_unit)

    def _get_matching_unit(self, matchobj):

        raw_unit = self[matchobj.group(2)]

        unit = raw_unit

        return unit

    def __setitem__(self, key, value):

        if isinstance(value, _pint_qty_type):
            unt = self._get_builtin_units(key)
            ref_value = value.to(unt).magnitude
        else:
            ref_value = value

        target_type = self._get_datatype(key)
        try:
            assert isinstance(ref_value, target_type)
        except AssertionError:
            ref_value = _type_recovery(ref_value, target_type)

        self._p_odr.utils.run_command('edit ' + self.fullname + ' ' + key + '=' + _odssrep(ref_value))

    def __str__(self):
        return '<PackedOpendssElement: {0}>'.format(self._name)

    def __repr__(self):
        return self.__str__()


class _CallFinalizer(Callable):
    def __init__(self, interface, selector_fns: tuple, name: str, s_interface_name: str):
        self._super_interface_name = s_interface_name
        self._interface = interface
        self._selectors = selector_fns
        self._name_to_select = name

    def __getattr__(self, item):
        return _CallFinalizer(getattr(self._interface, item), self._selectors, self._name_to_select,
                              self._super_interface_name)

    # no cache on call, otherwise it would not see recent edits
    def __call__(self, *args):
        # this is the key bit: the PackedOpendssElement that instances this class is capable of retaining its name and
        # auto select itself before calling the underlying odr.
        for sel in self._selectors:
            sel(self._name_to_select)

        logger.debug('Calling {0} with arguments {1}'.format(str(self._interface), str(args)))
        return self._interface(*args)

    @property
    def super_interface(self):
        return self._super_interface_name
