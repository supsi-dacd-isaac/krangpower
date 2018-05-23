# OpendssdirectEnhancer by Federico Rosato
# a wrapper for opendssdirect.py by Dheepak Krishnamurthy and Maximilian J. Zangs

import copy
import json
import logging
import os
import platform
import re
from functools import reduce
from logging.handlers import RotatingFileHandler as _RotatingFileHandler
from math import sqrt
from operator import getitem, attrgetter
from typing import Callable

import numpy as np
import opendssdirect as odr
from pandas import DataFrame
# from profilehooks import profile

from krangpower import components as co
from krangpower.aux_fcn import lower as _lower
from krangpower.aux_fcn import pairwise as _pairwise
from krangpower.components import resolve_unit, _SnpMatrix, _pint_qty_type, _odssrep, _type_recovery, _global_log_level
from krangpower.components import um, _mlog, config, thisdir

__all__ = ['OpendssdirectEnhancer']

# the default entity parameter values are loaded in order to allow correct type casting and comparison with the default
# when calling _Packed's __getitem__ method
_default_entities = co.default_comp

_classmap = co.get_classmap()

# information about the opendssdirect interface is loaded in order to feed the _PackedOpendssObject metaclass
_interface_methods_path = os.path.join(thisdir, config.get('data_files', 'interfaces'))
with open(_interface_methods_path, 'r') as ifile:
    itf = json.load(ifile)
_interface_methods = {(k,): v for k, v in itf.items()}
_unit_measurement_path = os.path.join(thisdir,  config.get('data_files', 'measurement_units'))
_treatments_path = os.path.join(thisdir, config.get('data_files', 'treatments'))
_interf_selectors_path = os.path.join(thisdir, config.get('data_files', 'interface_selectors'))
with open(_interf_selectors_path, 'r') as ifile:
    itf_sel_names = json.load(ifile)


# instantiating the module command logger
def _create_command_logger(name):
    logformat = '%(asctime)s - %(message)s'
    cmd_logger = logging.getLogger(name)
    cmd_logger.setLevel(_global_log_level)
    logformatter = logging.Formatter(logformat)

    # filehandler
    try:
        if platform.system() == 'Windows':
            logpath = os.path.join(os.getenv('APPDATA'), config.get('log_file', 'log_folder'), config.get('log_file', 'commands_log_path'))
        elif platform.system() == 'Linux':
            logpath = os.path.join('/var/log', config.get('log_file', 'log_folder'), config.get('log_file', 'commands_log_path'))
        else:
            raise OSError('Could not find a valid log path.')
        if not os.path.exists(os.path.dirname(logpath)):
            os.makedirs(os.path.dirname(logpath))
        maxsize = config.getfloat('log_settings', 'max_log_size_mb')
        fh = _RotatingFileHandler(logpath, maxBytes=maxsize*1e6, backupCount=0)
        fh.setFormatter(logformatter)
        fh.setLevel(logging.DEBUG)
        cmd_logger.addHandler(fh)
    except PermissionError:
        # this is handled to the console stream
        cmd_logger.warning('Permission to write log file denied')

    return cmd_logger


# The choice of a module-wide command logger is due to the fact that opendssdirect remains unique even when several
# Enhancers are instantiated. So they better refer to a unique command logger. The distinction between instances is
# possible through Enhancer.id, settable as kwarg oe_id when instantiating it.
_default_name = 'OpenDSSEnhancer'
_clog = _create_command_logger(_default_name)


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


def _cast_dumbstring(string: str, data_type):
    """Casts, if possible, a raw string returned by the OpenDSS text interface to the type specified in data_type."""

    if data_type == str:
        return string
    if data_type in (int, float):
        return data_type(string)
    elif data_type == np.matrix:
        return _SnpMatrix(string
                          .replace(' |', ';')
                          .replace('|', ';')
                          .replace('[', '')
                          .replace(' ]', '')
                          .replace(']', '')
                          .replace(' ', ','))
    elif data_type == list:
        dp_str = re.sub('[\,|\ ]*(\]|"|\))', '', string)
        dp_str = re.sub('(\[|"|\()\ *', '', dp_str)
        items = dp_str.split(',')
        try:
            return [int(x) for x in items]
        except ValueError:
            try:
                return [float(x) for x in items]
            except ValueError:
                return items
    else:
        raise TypeError('Could not cast the DSS property string "{1}": type {0} unknown'.format(str(data_type), string))


# there's no XYCurves.AllNames() or similar, so we have to mock up one ourselves
def _xycurve_names():

    if odr.XYCurves.Count() == 0:
        return []

    xynames = []
    i = 1
    odr.XYCurves.First()
    while i != 0:
        xynames.append(odr.XYCurves.Name())
        i = odr.XYCurves.Next()

    return xynames
# </editor-fold>


class OpendssdirectEnhancer:
    """OpendssdirectEnhancer constitues an OpenDSSDirect.py wrapper and overhaul, and it is usable in and by itself as
    a substitute of the OpenDSSDirect.py module when interested in the additional functionality.

    OpendssdirectEnhancer is designed to behave exactly as the opendssdirect MODULE by wrapping all its calls,
    but adds some handy functionality:

    -   Items returned as a list of floats (e.g., <>.Circuit.Losses(), <>.Bus.Voltages()...)
        are returned as lists of complex numbers, matrices of [nterm x ncond], etc. as is appropriate. How the results
        returned by opendss are treated is decided in the configuration files of krangpower.
    -   Structured items such as opendssdirect.Circuit.SistemY() are returned as pandas.DataFrame for easier
        manipulation, export...
    -   Items come, where appropriate, as Quantities (from the pint package) with the appropriate measurement unit.
        This enables easy conversions and secures against miscalculations.
    -   OpendssdirectEnhancer supports easy and intuitive element exploration through bracket indicization with the
        names, either fully qualified or simple, of the circuit's elements in the way shown below:

            >>> myOE = OpendssdirectEnhancer()
            >>> myOE.utils.run_command("Clear")
            ''
            >>> myOE.utils.run_command("New object = circuit.myckt bus1=sourcebus basekv=11.0 pu=1.0 angle=0.0 phases=3")
            ''
            >>> myOE.utils.run_command("New load.myload bus1=sourcebus kw=10.0 kv=11.0 basefreq=50.0")
            ''
            >>> myOE.utils.run_command("solve")
            ''
            >>> myOE['load.myload'].CFactor() # taken from the Loads interface
            4.0
            >>> abs(myOE['myload'].Currents()[0,1]) # taken from the CktElement interface
            <Quantity(0.5964385133615372, 'ampere')>

    As in the example above, the user needs not to worry about selecting the item with MOE.Circuit.SetActiveElement('load.myload'),
    MOE.Laods.Name('load.myload'); the OpendssdirectEnhancer takes care of the selections automatically.

    The OpendssdirectEnhancer members are dynamically generated from those of OpenDSSDirect.py, so the api reference is
    the same and can be found at (https://nrel.github.io/OpenDSSDirect.py/index.html)
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

    # loads chains of functions through which to pass the rough outputs of opendss.
    with open(_treatments_path, 'r') as tfile:
        rtrt = json.load(tfile)
    trt = dict()
    for subdic_name, subdic in rtrt.items():
        nsd = {k: tuple([globals()[t] for t in v]) for k, v in subdic.items()}
        trt[subdic_name] = nsd

    # loads measurement units for the interface of components without self-referencing.
    # the components with self referencing, like lines and loadshapes, are taken care of at runtime.
    with open(_unit_measurement_path, 'r') as ufile:
        rumr = json.load(ufile)
    umr = dict()
    for subdic_name, subdic in rumr.items():
        nsd = {k: um.parse_units(v) for k, v in subdic.items()}
        umr[subdic_name] = nsd

    _names_up2date = False
    _cached_allnames = []

    def __init__(self, stack=None, oe_id=_default_name):
        if stack is None:
            self.stack = []
        else:
            self.stack = stack
        self.id = oe_id

    def __getattr__(self, item):
        """Directly wraps the calls to OpenDSSDirect.py."""
        if self._chain_hasattr(self.stack + [item]):
            return OpendssdirectEnhancer(self.stack + [item], self.id)
            # return getattr(odr, item)
        else:
            raise AttributeError('Could not find the attribute {0} in the available interfaces.'.format(item))

    def __getitem__(self, item):
        """Bracket indicization looks for an object with the name desired and returns a nice packed element whose
         attributes are the same present in the OpenDSSDirect.py api, but without ever worrying again about
         SetActiveElement() and the likes.
        """

        try:
            assert item.lower() in map(lambda name: name.lower(), self.get_all_names())
        except AssertionError:
            bare_names_dict = {name.lower().split('.', 1)[1]: name.lower() for name in self.get_all_names()}
            try:
                assert item.lower() in bare_names_dict.keys()
            except AssertionError:
                raise KeyError('Element {0} was not found in the circuit'.format(item))
            else:
                fullitem = bare_names_dict[item.lower()]
        else:
            fullitem = item.lower()

        return _PackedOpendssElement(*fullitem.split('.', 1), self)

    def __str__(self):
        return '<OE:({0})>'.format('.'.join(['opendssdirect'] + self.stack[1:]))

    def __repr__(self):
        return self.__str__()

    def __call__(self, *args):
        # the dictionary of unit measures is found. The cases of lines and loadshapes are treated aside, because their
        # unit dictionary references the object itself ('units' or 'useactual'), so the dependency must be resolved
        # dynamically
        if self.stack[0] == 'Lines':
            um_d = self.line_umd(self.line_um[odr.Lines.Units()])
        elif self.stack[0] == 'LoadShape':
            um_d = self.loadshape_umd(self.line_um[bool(odr.LoadShape.UseActual())])
        else:
            um_d = self.umr

        # the dictionary is walked to find the unit of the particular value requested with this call
        try:
            ums = reduce(getitem, self.stack, um_d)  # nested dict search
        except KeyError:
            ums = None

        # the trt dictionary is walked to find which function must the output be passed through
        try:
            trt = reduce(getitem, self.stack, self.trt)
        except KeyError:
            trt = tuple()  # empty tuple in order to shortcut the following iteration when no treatments needed

        # we retrieve the desired member from opendssdirect.py and call it with the arguments
        odrobj = self._chain_getattr()
        e_ordobj = odrobj(*args)

        # an explicit control must be carried out over what's returned when we are calling the text interface.
        # this is because many errors from the text interface are returned silently as regular strings, potentially
        # leading to dangerous errors.

        # the result is finally treated and assigned a unit.
        for t in trt:
            e_ordobj = t(e_ordobj)

        return _assign_unit(e_ordobj, ums)

    @classmethod
    def get_all_names(cls):
        """Gets the fully qualified names of the elements, plus buses, loadshapes and xycurves.
        It's worth noting that object such as LineCodes, WireCodes, etc are not as of today retrievable, because their
        names are not accessible in a direct way."""
        odr.utils.run_command('makebuslist')
        if cls._names_up2date:
            return cls._cached_allnames
        else:
            anl = []
            odr.utils.run_command('makebuslist')
            anl.extend(map(lambda bn: 'bus.' + bn, odr.Circuit.AllBusNames()))
            anl.extend(odr.Circuit.AllElementNames())
            anl.extend(map(lambda ln: 'loadshape.' + ln, odr.LoadShape.AllNames()))
            anl.extend(map(lambda ln: 'xycurve.' + ln, _xycurve_names()))

            cls._names_up2date = True
            cls._cached_allnames = anl

            return anl

    @staticmethod
    def loadshape_umd(use_actual: bool):
        """Dynamically generates the measurement units dictionary for a loadshape based on the property use_actual."""

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

    def txt_command(self, cmd_str: str, echo=True):
        """Performs a text interface call with the argument passed and logs command and response. The log output can be
        suppressed by setting the keyword argument echo=False."""
        rslt = self.utils.run_command(cmd_str)  # rslt could be an error string too
        if echo:
            self.log_line('[' + cmd_str.replace('\n', '\n' + ' ' * (30 + len(_default_name)))
                          + ']-->[' + rslt.replace('\n', '') + ']')

        if self._influences_names(cmd_str):
            OpendssdirectEnhancer._names_up2date = False

        try:
            _validate_text_interface_result(rslt)
        except OpenDSSTextError:
            raise
        else:
            return rslt

    @staticmethod
    def _influences_names(cmd_str):
        if cmd_str.lower().startswith('new'):
            return True
        else:
            return False

    def log_line(self, line: str, lvl=logging.DEBUG):
        """Logs a line in the command log."""
        _clog.log(lvl, '(id:{0})-'.format(self.id) + line)

    def _chain_getattr(self, stack=None, start_obj=odr):
        if stack is None:
            stack = self.stack
        return attrgetter('.'.join(stack))(start_obj)

    def _chain_hasattr(self, stack=None, start_obj=odr):
        if stack is None:
            stack = self.stack

        if len(stack) == 1:
            return hasattr(start_obj, stack[0])
        else:
            return self._chain_hasattr(stack[1:], getattr(start_obj, stack[0]))


class OpenDSSTextError(Exception):
    """Meant to be thrown when the string returned by opendss text interface represents an error."""
    pass


def _validate_text_interface_result(result_string: str):
    """This function is passed the raw, direct output opendss text interface and performs checks on the results to see
    if the string returned is valid (and not, for example, a warning). This function either returns nothing or
    raises an error."""
    if result_string.lower().startswith(('warning', 'error')):
        raise OpenDSSTextError(result_string)

    if result_string.lower().find('not found') != -1:
        raise OpenDSSTextError(result_string)

    # this is what odr returns, instead of raising, if you request invalid props with ?
    if result_string == 'Property Unknown':
        raise KeyError('Property Unknown')


class _PackedOpendssElement:
    def __init__(self, eltype, name, p_odr):

        # reconstructs what are the available interfaces from file
        self._available_interfaces = tuple(getattr(p_odr, itf_name) for itf_name in itf_sel_names['interfaces'][eltype])

        # reconstructs from file what are the interfaces in which I have to select the element in order to get the
        # results that pertain it from self._available_interfaces
        self._selectors = []
        sel_attrs_chain = tuple(itf_sel_names['selectors'][eltype])
        for ac in sel_attrs_chain:
            for i, a in enumerate(ac.split('.')):
                if i == 0:
                    obi = getattr(p_odr, a)
                else:
                    obi = getattr(obi, a)
            self._selectors.append(obi)

        self._name = name
        # todo perform sanity checks on name
        # todo distinct linecode_a, etc and make sure that _eltype is coherent with components
        self._eltype = eltype
        self._p_odr = p_odr

        # dynamically add methods to expose
        for i in self._available_interfaces:
            for m in _interface_methods.get(tuple(i.__name__.split('.')[-1]), []):
                setattr(self, m, self._craft_member(m))

    @property
    def topological(self):
        """Returns those properties that are marked as 'topological' in the configuration files and identify the wiring
        location of the element. (Typically, bus1 and, if it exists, bus2.)"""
        top_par_names = _default_entities['default_' + self._eltype]['topological']

        rt = [self[t] for t in top_par_names.keys()]
        if len(rt) == 1 and isinstance(rt[0], list):
            return tuple(rt[0])
        else:
            return tuple(rt)

    @property
    def type(self):
        return self._eltype

    @property
    def fullname(self):
        return self._eltype + '.' + self._name

    @property
    def name(self):
        return self._name

    def _craft_member(self, item: str):
        """This second order function returns a function that invokes self._side_getattr on item. Such returned
         function are intended to be assigned to a _PackedOpendssElement's attribute with the same name as item, and
         such operation is carried out in __init__."""
        
        def pckdss_mthd(*args):
            # self._side_getattr(m) is a CallFinalizer
            return self._side_getattr(item)(*args)

        # poo_mthd.__repr__ = lambda: '<function {0} of {1}.{2}>'.format(item, self._eltype, self._name)
        # poo_mthd.__name__ = lambda: item
        return pckdss_mthd

    def _side_getattr(self, item):
        """Returns a _CallFinalizer pointing to item, if item is a valid interface identificator of the object in
        OpenDSSdirect.py.
        For example, If the object is a Isource, AngleDeg will be a valid item."""
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
        by calls to __getitem__."""
        try:
            props = self._p_odr[self.fullname].AllPropertyNames()
        except AttributeError:
            raise ValueError('{0}-type objects are not dumpable.'.format(self._eltype.upper()))
        return {p: self._p_odr[self.fullname][p] for p in props}

    # @profile(immediate=True)
    def unpack(self, verbose=False):
        """Returns a _DssEntity (or a descendant) corresponding to _PackedOpendssElement."""

        # identify the corresponding class in the components file
        myclass = _classmap[self._eltype]

        # properties are dumped
        all_props = self.dump()

        # the names of those properties that are ok to pass to the _DssEntity are taken from the components'
        # configuration file
        valid_props = copy.deepcopy(_default_entities['default_' + self._eltype]['properties'])
        valid_props.update(_default_entities['default_' + self._eltype].get('associated', {}))

        # either those properties dumped that are valid, or those that are valid AND different from the default values,
        # are stored in dep_prop
        if verbose:
            dep_prop = {k.lower(): v for k, v in all_props.items() if k.lower() in valid_props.keys()}
        else:
            dep_prop = {k.lower(): v for k, v in all_props.items() if k.lower() in valid_props.keys() and np.matrix(v != valid_props[k.lower()]).any()}

        # the _DssEntity is instantiated with the properties in dep_prop.
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

        # the raw property is returned in string form from Opendssdirect's text interface with the ? operator
        try:
            rslt = self._p_odr.utils.run_command('? {0}.{1}'.format(self.fullname, item))
        except KeyError:
            # the rudimental Keyerror retrieved by _validate... is further explained
            raise KeyError('Invalid property {1} requested for {0}'.format(self.fullname, item))

        # correct type casting
        try:
            data_type = self._get_datatype(item)
        except KeyError:
            # happens, for example, when topological parameters are requested, because while being valid they are not
            # listed in the default-type configuration files
            return rslt

        rslt = _cast_dumbstring(rslt, data_type)

        # an unit is possibly added
        unt = self._get_builtin_units(item)
        if unt is None:
            return rslt
        else:
            return rslt * unt

    def __pos__(self):
        """The + operator enables the _PackedOpendssElement, undoing the effects of a previous __neg__ call on the same
         element. It directly wraps the 'enable' command from opendss."""
        self._p_odr.txt_command('enable {0}'.format(self.fullname))

    def __neg__(self):
        """The - operator disables the _PackedOpendssElement, so that it's like it does not exists. It directly wraps the
        'disable' command from opendss."""
        self._p_odr.txt_command('disable {0}'.format(self.fullname))

    def _get_datatype(self, item):
        """Gets what type of data corresponds to the property name passed as argument. The information is retrieved from
        the configuration files."""
        try:
            return type(_default_entities['default_' + self._eltype]['properties'][item.lower()])
        except KeyError:
            return type(_default_entities['default_' + self._eltype]['topological'][item.lower()])

    def _get_builtin_units(self, item):
        """Gets what measurement unit corresponds to the property name passed as argument. The information is retrieved
        from the configuration files."""
        raw_unit = _default_entities['default_' + self._eltype]['units'].get(item.lower(), None)
        if raw_unit is None:
            return None
        else:
            return resolve_unit(raw_unit, self._get_matching_unit)

    def _get_matching_unit(self, matchobj):
        """Gets the property stored in the argument, a matchobj, group(2). It is a supporting function to
        _get_builtin_units, meant to retrieve, the property 'units' of the object or something similar."""
        raw_unit = self[matchobj.group(2)]
        return raw_unit

    def __setitem__(self, key, value):
        """Edits a property of the _PackedOpendssElement.
        Beware: If you do not pass a value with a pint measurement unit when editing a parameter that has a physical
        dimensionality, it will be assumed that you are using the default units."""

        # it is checked whether you passed a _pint_qty_type as value or not. Throughout the function, errors will be
        # trhown if: _pint_qty_type is passed for a property without unit, _pint_qty_type has the wrong dimensionality,
        # the content of the _pint_qty_type is not the right data type (such as a matrix instead of an int).
        if isinstance(value, _pint_qty_type):
            unt = self._get_builtin_units(key)
            ref_value = value.to(unt).magnitude
        else:
            ref_value = value

        # the target datatype is taken from the configuration files and checked
        target_type = self._get_datatype(key)
        try:
            assert isinstance(ref_value, target_type)
        except AssertionError:
            ref_value = _type_recovery(ref_value, target_type)

        # the edit is performed through the text interface with the 'edit' command
        self._p_odr.utils.run_command('edit ' + self.fullname + ' ' + key + '=' + _odssrep(ref_value))

    def __str__(self):
        return '<PackedOpendssElement: {0}>'.format(self._name)

    def __repr__(self):
        return self.__str__()


class _CallFinalizer(Callable):
    def __init__(self, interface, selector_fns: list, name: str, s_interface_name: str):
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

        _mlog.debug('Calling {0} with arguments {1}'.format(str(self._interface), str(args)))
        return self._interface(*args)

    @property
    def super_interface(self):
        return self._super_interface_name
