# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

# OpendssdirectEnhancer by Federico Rosato
# a wrapper for opendssdirect.py by Dheepak Krishnamurthy and Maximilian J. Zangs

import logging
import types
from copy import deepcopy as _deepcopy
from functools import reduce as _reduce, partial as _partial
from inspect import getmembers as _inspect_getmembers
from json import load as _json_load
from logging import DEBUG as _DBG_LVL
from math import sqrt as _sqrt
from operator import getitem as _getitem, attrgetter as _attrgetter
from re import sub as _sub
from sys import modules as _sys_modules

import numpy as _np
import opendssdirect as _odr
from pandas import DataFrame as _DataFrame

from .._aux_fcn import lower as _lower
from .._aux_fcn import pairwise as _pairwise
from .._components import LineGeometry
from .._components import _resolve_unit, _type_recovery, _odssrep, SnpMatrix
from .._components import get_classmap as _get_classmap
from .._config_loader import DEFAULT_ENH_NAME, UNIT_MEASUREMENT_PATH, TREATMENTS_PATH, ERROR_STRINGS, \
    UM as _UM, INTERFACE_METHODS_PATH, DEFAULT_COMP as _DEFAULT_COMP, PINT_QTY_TYPE, INTERF_SELECTORS_PATH
from .._logging_init import clog, mlog, get_log_level

try:
    assert _odr.Basic.Start(0)
except TypeError:
    # retrocompatibility with OpenDSSDirect.py <0.3
    assert _odr.Basic.Start()

# <editor-fold desc="Auxiliary functions">


def _assign_unit(item, unit: type(_UM.m) or None):
    if unit is None:
        return item
    elif isinstance(item, dict):
        return {k: v * unit for k, v in item.items()}
    elif isinstance(item, _DataFrame):
        # pandas' _DataFrame is a mess together with pint
        return item
    # elif hasattr(item, '__iter__'):
    #     # return [el * unit for el in item]
    #     return _asarray(item) * unit
    else:
        return item * unit


def _asarray(item):
    return _np.asarray(item)


def _couplearray(item):
    return _np.array(item[0::2]) + _np.array(item[1::2]) * 1j


def _terminalize_cpx(item):
    # note that, when I pass an item to terminalize, I am taking for granted that I can find nterm and ncond in the
    # respective calls to odr. If you called odr.CktElement.Powers(), for example, I take it that you knew what
    # you were doing. Calls coming from PackedElements, instead, should perform the cktelement selection just before
    # the call.

    nterm = _odr.CktElement.NumTerminals()
    ncond = _odr.CktElement.NumConductors()

    assert len(item) == nterm * ncond * 2
    cpxr = _np.zeros([nterm, ncond], dtype=_np.complex)

    for idx, couple in enumerate(_pairwise(item)):
        real = couple[0]
        imag = couple[1]
        cpxr[int(idx / ncond), (idx % ncond)] = _np.sum([_np.multiply(1j, imag), real], axis=0)

    return cpxr


def _terminalize_int(item):
    nterm = _odr.CktElement.NumTerminals()
    ncond = _odr.CktElement.NumConductors()

    assert len(item) == nterm * ncond * 1
    int_r = _np.zeros([nterm, ncond], dtype=_np.int)
    for idx, element in enumerate(item):
        int_r[int(idx / ncond), (idx % ncond)] = element

    return int_r


def _cpx(item):
    return item[0] + 1j*item[1]


def _dictionize_cktbus(item):
    return dict(zip(_lower(_odr.Circuit.AllBusNames()), item))


def _dictionize_cktels(item):
    return dict(zip(_lower(_odr.Circuit.AllElementNames()), item))


def _dictionize_cktnodes(item):
    return dict(zip(_lower(_odr.Circuit.YNodeOrder()), item))


def _matricize_ybus(item):

    raw_n_ord = _lower(_odr.Circuit.YNodeOrder())
    mtx = _matricize(item)
    return _DataFrame(data=mtx, index=raw_n_ord, columns=raw_n_ord)


def _matricize(item):

    if len(item) == 1:
        return item

    side = _sqrt(len(item)/2)
    assert side == int(side)
    side = int(side)
    mtx = _np.reshape(item[0::2], (side, side)) + \
          _np.reshape(item[1::2], (side, side)) * 1j

    return _np.transpose(mtx)  # it's transposed because it's originally given in column order


# there's no XYCurves.AllNames() or similar, so we have to mock up one ourselves
def _xycurve_names():

    if _odr.XYCurves.Count() == 0:
        return []

    xynames = []
    i = 1
    _odr.XYCurves.First()
    while i != 0:
        xynames.append(_odr.XYCurves.Name())
        i = _odr.XYCurves.Next()

    return xynames
# </editor-fold>


# <editor-fold desc="Loads and declarations">
_this_module = _sys_modules[__name__]
_classmap = _get_classmap()

_ID = 'OpendssdirectEnhancer'
setattr(_this_module, 'utils', _odr.utils)


# loads chains of functions through which to pass the rough outputs of opendss.
with open(TREATMENTS_PATH, 'r') as _tfile:
    _rtrt = _json_load(_tfile)
_trt = dict()
for _subdic_name, _subdic in _rtrt.items():
    _nsd = {k: tuple([globals()[_t] for _t in _v]) for k, _v in _subdic.items()}
    _trt[_subdic_name] = _nsd

# loads measurement units for the interface of components without self-referencing.
# the components with self referencing, like lines and loadshapes, are taken care of at runtime.
with open(UNIT_MEASUREMENT_PATH, 'r') as _ufile:
    _rumr = _json_load(_ufile)
_umr = dict()
for _subdic_name, _subdic in _rumr.items():
    _nsd = {_k: _UM.parse_units(_v) for _k, _v in _subdic.items()}
    _umr[_subdic_name] = _nsd


with open(INTERFACE_METHODS_PATH, 'r') as _ifile:
    _itf = _json_load(_ifile)
_interface_methods = {(_k,): _v for _k, _v in _itf.items()}
with open(INTERF_SELECTORS_PATH, 'r') as _ifile:
    _itf_sel_names = _json_load(_ifile)
# </editor-fold>


# <editor-fold desc="Text output check and conversion">
class OpenDSSTextError(Exception):
    """Meant to be thrown when the string returned by opendss text interface represents an error."""
    pass


# this ctxman is necessary for suppressing the rogue warnings from OpenDSSDirect.py
# when the situation is under control and they are unnecessary
class _LogKiller:
    def __enter__(self):
        logging.disable(logging.CRITICAL)

    def __exit__(self, a, b, c):
        logging.disable(logging.NOTSET)


def _validate_text_interface_result(result_string: str):
    """This function is passed the raw, direct output opendss text interface and performs checks on the results to see
    if the string returned is valid (and not, for example, a warning). This function either returns nothing or
    raises an error."""

    error_strings_begin = tuple(ERROR_STRINGS['beginning'])
    if result_string.lower().startswith(error_strings_begin):
        raise OpenDSSTextError(result_string)

    error_strings_middle = tuple(ERROR_STRINGS['middle'])
    for inner_message in error_strings_middle:
        if result_string.lower().find(inner_message) != -1:
            raise OpenDSSTextError(result_string)

    # this is what odr returns, instead of raising, if you request invalid props with ?
    if result_string == 'Property Unknown':
        raise KeyError('Property Unknown')


def _cast_dumbstring(string: str, data_type):
    """Casts, if possible, a raw string returned by the OpenDSS text interface to the type specified in data_type."""

    if data_type == str:
        return string
    if data_type in (int, float):
        if string == '----':  # this happens, f.e., if you define a line by rmatrix, xmatrix and then ask for r0, x0...
            return _np.NaN
        else:
            return data_type(string)
    elif data_type == _np.matrix:
        return SnpMatrix(string
                         .replace(' |', ';')
                         .replace('|', ';')
                         .replace('[', '')
                         .replace(' ]', '')
                         .replace(']', '')
                         .replace(' ', ','))
    elif data_type == list:
        dp_str = _sub('[\,|\ ]*(\]|"|\))', '', string)
        dp_str = _sub('(\[|"|\()\ *', '', dp_str)
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


def _influences_names(cmd_str):
    if cmd_str.lower().startswith('new'):
        return True
    else:
        return False
# </editor-fold>


# <editor-fold desc="Dynamic unitmeasure registries">
_line_um = {
    0: _UM.unitlength,
    1: _UM.mile,
    2: _UM.kft,
    3: _UM.km,
    4: _UM.m,
    5: _UM.ft,
    6: _UM.inch,
    7: _UM.cm
}


def _loadshape_umd(use_actual: bool):
    """Dynamically generates the measurement units dictionary for a loadshape based on the property use_actual."""

    if use_actual:
        return {'Loadshapes':
                    {'HrInterval': _UM.hr,
                     'MinInterval': _UM.min,
                     'PBase': _UM.kW,
                     'PMult': _UM.kW,
                     'QBase': _UM.kVA,
                     'QMult': _UM.kVA,
                     'SInterval': _UM.s,
                     'TimeArray': _UM.hr}}

    else:
        return {'Loadshapes':
                    {'HrInterval': _UM.hr,
                     'MinInterval': _UM.min,
                     'PBase': _UM.kW,
                     'PMult': _UM.dimensionless,
                     'QBase': _UM.kW,
                     'QMult': _UM.dimensionless,
                     'SInterval': _UM.s,
                     'TimeArray': _UM.hr}}


def _line_umd(unit_length):
    """Special case of call if the object is a line; needed because it's a type of object for which a "units"
    properties exists that influences the other quantities dimensions."""

    line_qties = {'Lines': {'C0': _UM.nF / unit_length,
                            'C1': _UM.nF / unit_length,
                            'CMatrix': _UM.nF / unit_length,
                            'EmergAmps': _UM.A,
                            'Length': unit_length,
                            'NormAmps': _UM.A,
                            'R0': _UM.ohm / unit_length,
                            'R1': _UM.ohm / unit_length,
                            'RMatrix': _UM.ohm / unit_length,
                            'Rg': _UM.ohm / unit_length,
                            'Rho': _UM.ohm / unit_length,
                            'X0': _UM.ohm / unit_length,
                            'X1': _UM.ohm / unit_length,
                            'XMatrix': _UM.ohm / unit_length,
                            'Xg': _UM.ohm / unit_length,
                            'Yprim': _UM.siemens / unit_length}}

    return line_qties
# </editor-fold>


# <editor-fold desc="Module cache variables">
names_up2date = False
_cached_allnames = []
# </editor-fold>


# <editor-fold desc="Dynamical module population">
def _enh_call(*args, stack, odrobj):
    # the dictionary of unit measures is found. The cases of lines and loadshapes are treated aside, because their
    # unit dictionary references the object itself ('units' or 'useactual'), so the dependency must be resolved
    # dynamically
    if stack[0] == 'Lines':
        um_d = _line_umd(_line_um[_odr.Lines.Units()])
    elif stack[0] == 'LoadShape':
        um_d = _loadshape_umd(_line_um[bool(_odr.LoadShape.UseActual())])
    else:
        um_d = _umr

    # the dictionary is walked to find the unit of the particular value requested with this call
    try:
        ums = _reduce(_getitem, stack, um_d)  # nested dict search
    except KeyError:
        ums = None

    # the trt dictionary is walked to find which function must the output be passed through
    try:
        ths_trt = _reduce(_getitem, stack, _this_module._trt)
    except KeyError:
        ths_trt = tuple()  # empty tuple in order to shortcut the following iteration when no treatments needed

    # we retrieve the desired member from opendssdirect.py and call it with the arguments
    e_ordobj = odrobj(*args)

    # an explicit control must be carried out over what's returned when we are calling the text interface.
    # this is because many errors from the text interface are returned silently as regular strings, potentially
    # leading to dangerous errors.

    # the result is finally treated and assigned a unit.
    for _t in ths_trt:
        e_ordobj = _t(e_ordobj)

    return _assign_unit(e_ordobj, ums)


class _FnWrp:
    def __init__(self, partial_fkt):
        self._undfcn = partial_fkt
        self._stk = '.'.join(partial_fkt.keywords['stack'])
        self.__name__ = self._stk
        self.name = self._stk

    def __call__(self, *args):
        return self._undfcn(*args)

    def __str__(self):
        return '<function OpendssdirectEnhancer.' + self._stk + '>'

    def __repr__(self):
        return self.__str__()


# dynamical population of the module to mock OpenDSSDirect
for _i1 in _inspect_getmembers(_odr.dss):
    _i1_name = _i1[0]
    if _i1_name in _itf.keys():
        setattr(_this_module, _i1_name, types.ModuleType(_this_module.__name__ + '.' + _i1_name))
    else:
        continue
    for _i2 in _inspect_getmembers(_i1[1]):
        if _i2[0].startswith('_'):  # this vulgar hack is to avoid special methods
            continue
        try:
            _stack = [_i1_name, _i2[0]]
            _frozencall = _partial(_enh_call, stack=_stack, odrobj=_attrgetter('.'.join(_stack))(_odr.dss))
        except (AttributeError,):
            continue

        try:
            setattr(getattr(_this_module, _i1_name), _i2[0], _FnWrp(_frozencall))
        except (TypeError, AttributeError):
            continue
# </editor-fold>


# mocks a standard selector, but for named entities
def _named_selector(name, eltype):
    nm = (eltype.lower(), name.lower())
    _this_module.Circuit.SetActiveClass(nm[0])
    _this_module.ActiveClass.Name(nm[1])
    assert _this_module.Element.Name().lower() == '.'.join(nm)


class _PackedOpendssElement:

    def __init__(self, eltype, name):

        # reconstructs what are the available interfaces from file
        self._available_interfaces = tuple(
            getattr(_this_module, itf_name) for itf_name in _itf_sel_names['interfaces'][eltype])

        # reconstructs from file what are the interfaces in which I have to select the element in order to get the
        # results that pertain it from self._available_interfaces
        self._selectors = []
        sel_attrs_chain = tuple(_itf_sel_names['selectors'][eltype])
        for ac in sel_attrs_chain:
            obi = _attrgetter(ac)(_this_module)
            self._selectors.append(obi)

        self._name = name
        # todo perform sanity checks on name
        self.eltype = eltype

        # dynamically add methods to expose
        for i in self._available_interfaces:
            # for m in _interface_methods.get(tuple(i._undobj.split('.')[-1]), []):
            for m in _interface_methods.get(tuple(i.__name__.split('.')[-1]), []):
                setattr(self, m, self._craft_member(m))

    @property
    def topological(self):
        """Returns those properties that are marked as 'topological' in the configuration files and identify the wiring
        location of the element. (Typically, bus1 and, if it exists, bus2.)"""
        try:
            top_par_names = _DEFAULT_COMP['default_' + self.eltype]['topological']
        except KeyError:
            return None

        rt = [self[t] for t in top_par_names.keys()]
        if len(rt) == 1 and isinstance(rt[0], list):
            return tuple(rt[0])
        else:
            return tuple(rt)

    @property
    def type(self):
        return self.eltype

    @property
    def fullname(self):
        return self.eltype + '.' + self._name

    @property
    def name(self):
        return self._name

    @property
    def help(self):
        """Displays informations about the object's parameters."""
        return self.unpack().paramhelp

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
                return _CallFinalizer(getattr(itf, item), self._selectors, self.fullname, str(itf))
                # break unnecessary
            else:
                continue
        else:
            raise AttributeError('"{0}" is neither an attribute of the _PackedOpendssElement nor of the {1}-type object'
                                 ' it wraps.'
                                 .format(item, self.eltype.upper()))

    def dump(self):
        """Returns a dict with all the properties-values pairs of the object as they would be returned
        by calls to __getitem__."""
        try:
            for sel in self._selectors:
                try:
                    sel(self.fullname)
                except TypeError:
                    sel(self.name, self.eltype)  # named selector
            props = _this_module.Element.AllPropertyNames()
        except AttributeError:
            raise ValueError('{0}-type objects are not dumpable.'.format(self.eltype.upper()))
        return {p: self[p] for p in props}
    # @profile(immediate=True)

    def unpack(self, verbose=False):
        """Returns a _DssEntity (or a descendant) corresponding to _PackedOpendssElement."""

        # identify the corresponding class in the components file

        # the linegeometry internal structure is peculiar and has to be treated differently
        if self.eltype == 'linegeometry':
            return _unpack_linegeom(self)

        myclass = _classmap[self.eltype]

        # properties are dumped
        all_props = self.dump()

        # the names of those properties that are ok to pass to the _DssEntity are taken from the components'
        # configuration file
        valid_props = _deepcopy(_DEFAULT_COMP['default_' + self.eltype]['properties'])
        valid_props.update(_DEFAULT_COMP['default_' + self.eltype].get('associated', {}))

        ignored_props = _DEFAULT_COMP['default_' + self.eltype].get('ignored', [])
        redundant_props = _DEFAULT_COMP['default_' + self.eltype].get('redundant', [])

        valid_props = {k: v for k, v in valid_props.items() if k not in ignored_props + redundant_props}

        # either those properties dumped that are valid, or those that are valid AND different from the default values,
        # are stored in dep_prop
        if verbose:
            dep_prop = {k.lower(): v for k, v in all_props.items() if k.lower() in valid_props.keys()}
        else:
            dep_prop = {k.lower(): v for k, v in all_props.items() if
                        k.lower() in valid_props.keys() and _np.matrix(v != valid_props[k.lower()]).any()}

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
            rslt = txt_command('? {0}.{1}'.format(self.fullname, item))
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
        _this_module.txt_command('enable {0}'.format(self.fullname))

    def __neg__(self):
        """The - operator disables the _PackedOpendssElement, so that it's like it does not exists. It directly wraps the
        'disable' command from opendss."""
        _this_module.txt_command('disable {0}'.format(self.fullname))

    def _get_datatype(self, item):
        """Gets what type of data corresponds to the property name passed as argument. The information is retrieved from
        the configuration files."""
        try:
            return type(_DEFAULT_COMP['default_' + self.eltype]['properties'][item.lower()])
        except KeyError:
            return type(_DEFAULT_COMP['default_' + self.eltype]['topological'][item.lower()])

    def _get_builtin_units(self, item):
        """Gets what measurement unit corresponds to the property name passed as argument. The information is retrieved
        from the configuration files."""
        raw_unit = _DEFAULT_COMP['default_' + self.eltype]['units'].get(item.lower(), None)
        if raw_unit is None:
            return None
        else:
            return _resolve_unit(raw_unit, self._get_matching_unit)

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
        # thrown if: _pint_qty_type is passed for a property without unit, _pint_qty_type has the wrong dimensionality,
        # the content of the _pint_qty_type is not the right data type (such as a matrix instead of an int).
        if isinstance(value, PINT_QTY_TYPE):
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
        txt_command('edit ' + self.fullname + ' ' + key + '=' + _odssrep(ref_value))

    def __str__(self):
        return '<PackedOpendssElement({0})>'.format(self.fullname)

    def __repr__(self):
        return self.__str__()


class _CallFinalizer:

    def __init__(self, interface, selector_fns: list, name: str, s_interface_name: str):
        self._super_interface_name = s_interface_name
        self._interface = interface
        self._selectors = selector_fns
        self._name_to_select = name.split('.', 1)[1]
        self._eltype = name.split('.', 1)[0]  # it only is useful for the named_selector

    def __getattr__(self, item):
        return _CallFinalizer(getattr(self._interface, item), self._selectors, self._name_to_select,
                              self._super_interface_name)
    # no cache on call, otherwise it would not see recent edits

    def __call__(self, *args):
        # this is the key bit: the PackedOpendssElement that instances this class is capable of retaining its name and
        # auto select itself before calling the underlying odr.
        for sel in self._selectors:
            try:
                sel(self._name_to_select)
            except TypeError:
                sel(self._name_to_select, self._eltype)

        mlog.debug('Calling {0} with arguments {1}'.format(str(self._interface), str(args)))
        return self._interface(*args)

    @property
    def super_interface(self):
        return self._super_interface_name


def _unpack_linegeom(pckob):
    nc = pckob['nconds']

    nmc = {cc.split('.', 1)[1]: cc.split('.', 1)[0] for cc in get_all_names()}

    x = []
    h = []
    wrcn = []

    for cd in range(nc):
        pckob['cond'] = cd + 1
        cabtype = nmc[pckob['wire']]
        wrcn.append(pckob['wire'])
        x.append(pckob['x'].to('m').magnitude[0, 0])
        h.append(pckob['h'].to('m').magnitude[0, 0])

    naive_props = pckob.dump()
    del naive_props['x']
    del naive_props['h']
    del naive_props['wire']
    del naive_props['wires']
    del naive_props['cncable']
    del naive_props['cncables']
    del naive_props['tscable']
    del naive_props['tscables']
    del naive_props['units']
    del naive_props['normamps']
    del naive_props['emergamps']
    del naive_props['like']

    naive_props['x'] = _np.asmatrix(x)
    naive_props['h'] = _np.asmatrix(h)
    naive_props[cabtype] = wrcn

    return LineGeometry(pckob.name, **naive_props)


# <editor-fold desc="Exposed functions">
def pack(item):
    """Returns a PackedOpendssElement corresponding to item."""
    try:
        assert item.lower() in map(lambda name: name.lower(), _this_module.get_all_names())
    except AssertionError:
        bare_names_dict = {name.lower().split('.', 1)[1]: name.lower() for name in _this_module.get_all_names()}
        try:
            assert item.lower() in bare_names_dict.keys()
        except AssertionError:
            raise KeyError('Element {0} was not found in the circuit'.format(item))
        else:
            fullitem = bare_names_dict[item.lower()]
    else:
        fullitem = item.lower()

    return _PackedOpendssElement(*fullitem.split('.', 1))


def get_all_names():
    """Gets the fully qualified names of the elements, plus buses, loadshapes and xycurves."""
    _odr.utils.run_command('makebuslist')
    if _this_module.names_up2date:
        return _this_module._cached_allnames
    else:
        anl = []
        _odr.utils.run_command('makebuslist')
        anl.extend(map(lambda bn: 'bus.' + bn, _odr.Circuit.AllBusNames()))
        anl.extend(_odr.Circuit.AllElementNames())
        anl.extend(map(lambda ln: 'loadshape.' + ln, _odr.LoadShape.AllNames()))
        anl.extend(map(lambda ln: 'xycurve.' + ln, _xycurve_names()))

        # we access here all the elements that are interfaced through ActiveClass
        plus_classes = ('tsdata', 'cndata', 'linecode', 'wiredata', 'linegeometry')
        with _LogKiller():  # opendssdirect.py sends rogue warnings when an empty array is selected
            for pc in plus_classes:
                _this_module.Circuit.SetActiveClass(pc)
                anl.extend([pc + '.' + x for x in _this_module.ActiveClass.AllNames()])

        anl = [x.lower() for x in anl]

        _this_module.names_up2date = True
        _this_module._cached_allnames = anl

        return anl


def txt_command(cmd_str: str, echo=True):
    """Performs a text interface call with the argument passed and logs command and response. The results are checked for
     silent errors. **When instantiating components through this function, the update of the names returned by
     get_all_names()is triggered**. The log output can be suppressed by setting the keyword argument echo=False; but even
     in this case, if kp.get_log_level() is 0, the log of the command will be forced."""

    #with stdout_redirected2():
    rslt = _this_module.utils.run_command(cmd_str)  # rslt could be an error string too
    if echo or get_log_level() == 0:
        log_line('[' + cmd_str.replace('\n', '\n' + ' ' * (30 + len(DEFAULT_ENH_NAME)))
                 + ']-->[' + rslt.replace('\n', '') + ']')

    if _influences_names(cmd_str):
        _this_module.names_up2date = False

    try:
        _validate_text_interface_result(rslt)
    except OpenDSSTextError:
        raise
    else:
        return rslt


def log_line(line: str, lvl=_DBG_LVL):
    """Logs a line in the command log."""
    clog.log(lvl, '(id:{0})-'.format(_ID) + line)
# </editor-fold>


# <editor-fold desc="Rolling over the logs">
for lh in mlog.handlers:
    if hasattr(lh, 'doRollover'):
        lh.doRollover()

for lh in clog.handlers:
    if hasattr(lh, 'doRollover'):
        lh.doRollover()
# </editor-fold>
