import configparser as cfp
import copy
import csv
import json
from logging.handlers import RotatingFileHandler
import os.path
import re
import xml.etree.ElementTree as Xe
from abc import abstractmethod
from collections import OrderedDict, namedtuple
from functools import lru_cache
from sys import modules
import logging

import matplotlib.pyplot as plt
import numpy as np
import opendssdirect as odr
import scipy.io as sio
from dateutil.parser import parse as dateparse
from opendssdirect.utils import run_command as engine_command
from pandas import DataFrame, read_csv

import utils.aux_fcn
import pint


__all__ = ['load_entities', 'load_dictionary_json', 'CsvLoadshape', 'LineGeometry_C', 'LineGeometry_T', 'LineGeometry_O',
           'LineCode_A', 'LineCode_S', 'Line', 'WireData', 'CNData', 'TSData', 'Curve', 'PtCurve', 'EffCurve', 'Vsource',
           'Isource', 'DecisionModel', 'Load', 'Transformer', 'Capacitor', 'Capcontrol', 'Regcontrol','Reactor',
           'Monitor', 'BusVoltageMonitor', 'StorageController', 'Storage', 'PvSystem', 'FourQ', 'default_comp', 'logpath']

# from .odd import opendssdirect_treat
# odr = opendssdirect_treat()

um = pint.UnitRegistry()
um.define('percent = 0.01 * dimensionless = pct')
um.define('none = [generic_length] = unitlength')  # when lengths are set as none, this creates a common basis
um.define('mt = meter')

_pint_qty_type = type(1 * um.m)



# OPENDSS WRAPPER - WRITTEN BY FEDERICO ROSATO
# -------------------------------------------------------------

# -------------------------------------------------------------
# CONFIG LOAD
# -------------------------------------------------------------

thisdir = os.path.dirname(os.path.realpath(__file__))
config = cfp.ConfigParser()
config.read(os.path.join(thisdir, 'krang.cfg'))
default_settings_path_config = os.path.join(thisdir, config.get('data_files', 'default_settings'))
default_entities_path = os.path.join(thisdir, config.get('data_files', 'default_entities'))


# -------------------------------------------------------------
# MAIN LOGGER
# -------------------------------------------------------------

logformat = '%(asctime)s - %(levelname)s (%(funcName)s) -------------> %(message)s'
logger = logging.getLogger('krangpower')
logger.setLevel(logging.DEBUG)
logformatter = logging.Formatter(logformat)

# streamhandler
_ch = logging.StreamHandler()
_ch.setLevel(logging.WARN)
_ch.setFormatter(logformatter)
logger.addHandler(_ch)

# filehandler
try:
    logpath = os.path.join(os.getenv('APPDATA'), config.get('log_file', 'path'))
    if not os.path.exists(os.path.dirname(logpath)):
        os.makedirs(os.path.dirname(logpath))
    fh = RotatingFileHandler(logpath, maxBytes=2e6, backupCount=0)
    fh.setFormatter(logformatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
except PermissionError:
    # this is handled to the console stream
    logger.warning('Permission to write log file denied')



# <editor-fold desc="AUX FUNCTIONS">

def _cpx(odss_tuple, nterm, ncond):
    """
    This function transforms the raw data for electric parameters (voltage, current...) in a suitable complex array

    :param odss_tuple: tuple of nphases*2 floats (returned by odsswr as couples of real, imag components, for each phase of each terminal)
    :type odss_tuple: tuple or list
    :param nterm: number of terminals of the underlying electric object
    :type nterm: int
    :param ncond: number of conductors per terminal of the underlying electric object
    :type ncond: int
    :returns: a [nterm x ncond] numpy array of complex floats
    :rtype: numpy.ndarray
    """
    assert len(odss_tuple) == nterm * ncond * 2
    cpxr = np.zeros([nterm, ncond], 'complex')

    def pairwise(iterable):
        # "s -> (s0, s1), (s2, s3), (s4, s5), ..."
        a = iter(iterable)
        return zip(a, a)

    for idx, couple in enumerate(pairwise(odss_tuple)):
        real = couple[0]
        imag = couple[1]
        cpxr[int(idx / ncond), (idx % ncond)] = np.sum([np.multiply(1j, imag), real], axis=0)
        cpxr[int(idx / ncond), (idx % ncond)] = np.sum([np.multiply(1j, imag), real], axis=0)

    return cpxr


def _odssrep(data_raw):
    """
    This function, useful during the OpenDSS code generation, takes the stored parameters of the electrical components
    (such as impedance matrices, etc) as input and returns a string that represent that data according to the odsswr
    syntax
    :param data: the data to be converted
    :rtype: str
    """

    if isinstance(data_raw, _pint_qty_type):
        data = data_raw.magnitude
    else:
        data = data_raw

    def issym(matrix):
        return (np.rot90(np.flipud(np.tril(matrix)), -1) == np.triu(matrix)).all()

    def istriangular(number):
        #  remember that a number is triangular if T =n(n+1)/2 with n integer
        n = (-1 + np.sqrt(1 + 8 * number)) / 2
        return n == int(n)

    if isinstance(data, (float, int, str)):
        return str(data)

    elif isinstance(data, _SnpMatrix):

        order = data.diagonal().size
        ss = ''

        if issym(data):
            ind = np.tril_indices(order)  # returns tuple with 2 int np.arrays containing the indices

            for k, (r, c) in enumerate(zip(*ind)):
                ss += str(data[r, c]) + " "
                if istriangular(k + 1):  # after every triangular number, insert a return
                    ss += "| "
        else:
            for index, x in np.ndenumerate(data):
                ss += str(x) + " "
                if index[1] == order - 1:  # that means you're at the last item of that row
                    ss += "| "

        return '[' + ss[:-2] + ']'  # shaves last "| "

    elif isinstance(data, np.matrix) and data.shape[0] != 1:

        order = data.diagonal().size
        ss = ''

        if issym(data):
            ind = np.tril_indices(order)  # returns tuple with 2 int np.arrays containing the indices

            for k, (r, c) in enumerate(zip(*ind)):
                ss += str(data[r, c]) + " "
                if istriangular(k + 1):  # after every triangular number, insert a return
                    ss += "| "
        else:
            for index, x in np.ndenumerate(data):
                ss += str(x) + " "
                if index[1] == order - 1:  # that means you're at the last item of that row
                    ss += "| "

        return '[' + ss[:-2] + ']'  # shaves last "| "

    elif isinstance(data, np.matrix) and data.shape[0] == 1:  # i.e, a floatarray or intarray
        return str(data)[1:-1]

    elif isinstance(data, list):
        return '[' + ' '.join([str(d) for d in data]) + ']'

    else:
        raise TypeError


def _termrep(terminals):
    """
    This function takes a terminal collection (represented by a tuple of ints) and returns a representation that can be
    cat to a bus name in order to form a full bus-terminal qualification according to the odsswr syntax.
    Example: (1, 3, 2) ----> ".1.3.2"
    :param terminals: tuple of ints
    :type terminals: tuple
    :rtype: string
    """
    if terminals is None:
        return ''
    else:
        s = '.'
        try:
            for t in terminals:
                s += str(t) + '.'  # not very pythonic
            return s[0:-1]  # shaves final dot
        except TypeError:  # todo clean
            return '.' + str(terminals)


def _prcx(complex_number, prec=4, color=False):
    if color:
        return '\033[94m' + str(round(np.abs(complex_number), prec)) + '\033[91m' + '∠' + str(
            round(np.rad2deg(np.angle(complex_number)), prec)) + '°' + '\033[0m'
    else:
        return str(round(np.abs(complex_number), prec)) + '∠' + \
               str(round(np.rad2deg(np.angle(complex_number)), prec)) + '°'


def _is_timestamp(item):
    try:
        dateparse(item)
        return True
    except ValueError:
        return False


def _is_numeric_data(item):
    return re.fullmatch('([0-9]|\,|\.| )*', item) is not None


def _matricize(value):

    def desym(lol):
        size = len(lol)
        dsm = np.matrix(np.zeros([size, size]))
        for r in range(size):
            for c in range(r+1):
                dsm[r, c] = lol[r][c]
                dsm[c, r] = lol[r][c]

        return dsm

    if isinstance(value[0], str):
        return value
    else:
        try_mtx = np.matrix(value)
        if try_mtx.dtype == 'object':
            return desym(value)
        else:
            return try_mtx

# </editor-fold>


# <editor-fold desc="AUX CLASSES">
class _StringArray(list):
    """List intended to be initialized by a single string containing comma-separated values and to be able to __repr__
    that same string, albeit implementing list manipulation and indicization. Useful for manipulating and printing
    odsswr comma-separated values"""

    def __init__(self, datastring):
        super().__init__()
        for s in datastring.split(','):
            self.append(s)

    def __repr__(self):
        fstr = ''
        for s in self:
            fstr += s + ','
        return fstr[:-1]

    def __str__(self):
        return self.__repr__()


class _SnpMatrix(np.matrix):  # extends np.matrix, allowing to instantiate a symmetrical mtx by passing a tril string.
    """_SnpMatrix extends numpy.matrix. numpy.matrix can be initialized by a 'v11,v12;v21,v22'- like string; _SnpMatrix,
    in addition to this, can be initialized by a 'v11;v21,v22'-like string, representing the tril of a symmetrical
    matrix."""

    def __new__(subtype, data, dtype=None, copy=True):
        try:
            return super().__new__(subtype, data, dtype, copy)
        except ValueError as e:
            if dtype is None:
                c_dtype = float
            else:
                c_dtype = dtype

            if str(e) == 'Rows not the same size.' and isinstance(data, str):
                x = [c_dtype(l) for l in data.replace(',', ';').split(';')]
                order_r = (-1 + np.sqrt(1 + 8 * len(x))) / 2
                order = int(order_r)
                assert order == order_r  # if not, sthg is horribly wrong
                matrox = np.zeros([order, order])
                ind = np.tril_indices(order)  # returns tuple with 2 int np.arrays containing the indices
                for k, (r, c) in enumerate(zip(*ind)):
                    matrox[r, c] = x[k]
                    if r != c:  # don't write two times on the diagonal
                        matrox[c, r] = x[k]
                return super().__new__(subtype, matrox, c_dtype)
            else:
                raise
# </editor-fold>


ai_record = namedtuple('ai_record', ['element', 'bus', 'name'])


def load_dictionary_json(path):
    this_module = modules[__name__]
    classmap = {}
    for item in dir(this_module):
        classmap[item.lower()] = getattr(this_module, item)

    with open(path, 'r') as file:
        dik = json.load(file)

    # json entity file contain jsonized objects. This means that all lists are np.matrix.tolist representation
    # and we have to convert them back.
    for entity in dik:
        for prop, value in dik[entity]['properties'].items():
            if isinstance(value, list):
                dik[entity]['properties'][prop] = _matricize(value)

            # todo give it a unit measure

    return dik


def get_classmap():

    this_module = modules[__name__]
    classmap = {}
    for item in dir(this_module):
        classmap[item.lower()] = getattr(this_module, item)

    return classmap


def load_entities(path):

    classmap = get_classmap()

    with open(path, 'r') as file:
        dik = json.load(file)

    # json entity file contain jsonized objects. This means that all lists are np.matrix.tolist representation
    # and we have to convert them back.
    for entity in dik:
        for property, value in dik[entity]['properties'].items():
            if isinstance(value, list):
                dik[entity]['properties'][property] = _matricize(value)

    dicky = {}

    for entity_name in dik:
        elcls = classmap[dik[entity_name]['type']]
        if elcls.isnamed():
            dicky[entity_name] = elcls(entity_name, xml_rep=dik[entity_name]['properties'])
        else:
            dicky[entity_name] = elcls(dik[entity_name]['properties'])

    return dicky


def _resolve_unit(ustring: str, match_unit_getter):

    # resolution of self-referencing units, of the type "ohm / this.units" where "units" is a string parameter
    # of the _DSSObj
    if not isinstance(ustring, (list, tuple)):
        istring = [ustring]
    else:
        istring = ustring

    units = []

    for uidx, string in enumerate(istring):
        old_string = ''
        while old_string != string:
            old_string = string
            try:
                string = re.sub('(this\.)([^ ]+)', match_unit_getter, string)
            except TypeError:
                # unit getter returned a list
                def ginew(arg):
                    return match_unit_getter(arg, indx=uidx)
                string = re.sub('(this\.)([^ ]+)', ginew, string)

        units.append(string)

    if len(units) == 1:
        return [um.parse_units(s) for s in units][0]
    else:
        return [um.parse_units(s) for s in units]


def _type_recovery(value, target_type):
    try:
        if isinstance(value, int):
            assert target_type == float
            recovered_value = float(value)

        elif isinstance(value, list):
            assert target_type == np.matrix
            recovered_value = np.matrix(value)
        else:
            raise AssertionError
    except AssertionError:
        raise TypeError('Unrecoverable type {0}-->{1}'
                        .format(
                                type(value),
                                target_type))
    return recovered_value

# -------------------------------------------------------------
# GLOBAL DEFAULT ENTITIES DICT
# -------------------------------------------------------------
default_comp = load_dictionary_json(default_entities_path)


# -------------------------------------------------------------
# GENERIC DSSENTITY AND ITS NAMED VERSION
# -------------------------------------------------------------

class _DSSentity:  # implements the dictionary param, the xmlc drive for load and the on-the-fly overwrite
    
    _softmuldict = {'generator': {'CsvLoadshape': 'duty'},
                    'load': {'CsvLoadshape': 'duty'},
                    'vsource': {'CsvLoadshape': 'duty'},
                    'isource': {'CsvLoadshape': 'duty'},
                    'fourq': {'CsvLoadshape': 'duty',
                              'DecisionModel': '_dm'},
                    'linegeometry_c': {'CNData': 'cncable'},
                    'linegeometry_t': {'TSData': 'tscable'},
                    'linegeometry_o': {'WireData': 'wire'},
                    'line': {'LineGeometry_C': 'geometry',
                             'LineGeometry_O': 'geometry',
                             'LineGeometry_T': 'geometry',
                             'LineCode_S': 'linecode',
                             'LineCode_A': 'linecode'
                             },
                    'pvsystem': {'CsvLoadshape': 'duty',
                                 'PtCurve': 'p-tcurve',
                                 'EffCurve': 'effcurve'}
                    }

    _hardmuldict = {'Regcontrol': {'element': ('Transformer',)},
                    'Capcontrol': {'element': ('Capacitor',)},
                    'StorageController': {'element': ('CircuitElement',),
                                          'elementlist': ('list', 'Storage')},
                    'Monitor': {'element': ('_CircuitElementNBus',)}
                    }

    @classmethod
    def isnamed(cls):
        return False

    @classmethod
    def isai(cls):
        return False

    def __mul__(self, other):
        if self.toe in [x.lower() for x in self._softmuldict.keys()]:

            mytypemap = self._softmuldict[self.toe]

            if isinstance(other, list):
                datatype = other[0].__class__.__name__
                assert datatype in mytypemap.keys()
                assert all([x.__class__.__name__ == datatype for x in other])
                value = [x.name for x in other]
            else:
                datatype = other.__class__.__name__
                value = other.name

            prop_to_set = mytypemap[datatype]
            assert prop_to_set in self._associated.keys()
            self[prop_to_set] = value

            return self

        elif self.toe in [x.lower() for x in self._hardmuldict.keys()]:
            mypropmap = self._softmuldict[self.toe]
            for prop, proptype in mypropmap.items():

                g = other
                for chaintype in proptype:  # broken for storagecon
                    assert g.__class__.__name__ == chaintype
                    g = g[0]

                if isinstance(other, list):
                    value = [x.name for x in other]
                else:
                    value = other.name

                self[prop] = value
            return self

        else:
            raise TypeError

    def __init__(self, xml_rep, **kwargs):
        self.term_perm = None
        self.name = ''

        self.toe = self.__class__.__name__.lower()
        self.params = None
        self.editedParams = None
        self.default_units = None
        self._load_default_parameters()
        self.last_edited = []
        if not (xml_rep is None):
            try:
                self._setparameters(**xml_rep.load())
            except AttributeError:
                self._setparameters(**xml_rep)
        # todo write automation of typization (don't force the user to use types such as np.matrix, etc)
        self._setparameters(**kwargs)
        self.logger = logging.Logger('dummy')  # circuit's addc gives the component a logger after initializazion

    def __call__(self, **kwargs):
        # a call returns a copy of the object edited on the fly with the kwargs passed
        edited = copy.deepcopy(self)
        edited._setparameters(**kwargs)
        return edited

    # def __getattr__(self, item):
    #     # fallback getattr tries to call getparameters
    #     try:
    #         tp = self.params[item]
    #     except KeyError:
    #         raise AttributeError
    #
    #     return tp

    def edit(self, **kwargs):
        """
        Edits the parameters passed by kwargs. Requirements:

        - The parameter name specified is among the characteristic parameters of the object. See the documentation for
          the specific object for a review of the available parameters and their types.
        - The parameter value is of the appropriate type

        Passes when called with no arguments.
        """
        if kwargs:  # if arguments are provided, this is just an alias for setparameters
            self._setparameters(**kwargs)
        else:
            self.logger.info('Edit was invoked blank on a %s', self.toe)

        self.last_edited = list(kwargs.keys())

    def __getitem__(self, item):
        try:
            return self._getparameters(item)
        except KeyError:
            try:
                return self._associated[item]
            except KeyError:
                raise KeyError('Item {0} was not found in own parameters or associated entities'.format(item))

    def __setitem__(self, key, value):

        #move in setpar

        try:
            self._setparameters(**{key: value})
        except AttributeError:
            try:
                self._associated[key] = value
            except KeyError:
                raise KeyError('Item {0} was not found in own parameters or associated entities'.format(item))

    # def p(self, param, *default):  # just a getParameters alias for better usability
    #     """
    #     Gets the object parameter specified by param. A default return value can be specified in *default, should the
    #     original parameter be None. Throws KeyError if an invalid parameter is requested.
    #     DEPRECATED IN FAVOUR OF DEFAULT __getattr__
    #
    #     :param param: the parameter name
    #     :type param: string
    #     """
    #
    #     tp = self._getparameters(param)
    #     if tp is not None or len(default) == 0:
    #         return tp
    #     else:
    #         if len(default) > 1:
    #             raise IndexError
    #         return default[0]

    def _setparameters(self, **kwargs):

        def unitsfirst(dicky):
            for k, v in dicky.items():
                if k in ('units', 'runits', 'gmrunits'):
                    yield k, v
                else:
                    continue

            for k, v in dicky.items():
                if k not in ('units', 'runits', 'gmrunits'):
                    yield k, v
                else:
                    continue

        for parameter_raw, value_raw in unitsfirst(kwargs):

            # patch for the 'pct in keyword' problem
            if re.match('pct', parameter_raw.lower()):
                parameter = parameter_raw.lower().replace('pct', '%')
            else:
                parameter = parameter_raw

            if parameter.lower() in self.params.keys():
                target_list = self.params
                default_dict = 'properties'
            elif parameter.lower() in self._associated.keys():
                default_dict = 'associated'
                target_list = self._associated

            else:
                raise AttributeError('Tried to set unknown parameter {0}.'.format(parameter))
                # self.logger.warning('Tried to set unknown parameter %s. Blatantly ignored.', parameter)
                # should this raise an exception instead?

            # pint quantity check and conversion
            if isinstance(value_raw, _pint_qty_type):
                unt = _resolve_unit(self.default_units[parameter], self._get_matching_unit)
                if unt == um.none:
                    pass
                    # assert parameter_raw == 'length'
                    # unt = value_raw.units
                    # self['units'] = str(unt)  # todo this is wrong
                value = value_raw.to(unt).magnitude
            else:
                value = value_raw

            # type-correctness check of the raw value
            try:
                assert isinstance(value, (self.params_types_raw[parameter.lower()], type(None)))
            except AssertionError:
                target_type = self.params_types_raw[parameter.lower()]
                try:
                    value = _type_recovery(value, target_type)

                # if the value couldn't be salvaged, raise
                except TypeError:
                    raise TypeError('parameter "{0}" is of type {1} instead of {2} and could not be converted'
                                    .format(parameter.lower(),
                                            type(value),
                                            target_type))

            if isinstance(value, np.matrix):
                test = np.array_equal(value, default_comp['default_' + self.toe][default_dict][parameter])
            elif isinstance(value, list):
                test = value == default_comp['default_' + self.toe][default_dict][parameter]
            else:
                test = value == default_comp['default_' + self.toe][default_dict][parameter]

            if test:
                logger.debug('[{2}-{3}]Ignored setting {0} = {1} because identical to default'.format(parameter, str(value), self.toe, self.name))
                return

            # finally setting the parameter
            target_list[parameter.lower()] = value
            self.editedParams.append(parameter.lower())

    def _getparameters(self, param):

        if param.lower() in self.params.keys():
            target_list = self.params
        elif param.lower() in self._associated.keys():
            target_list = self._associated
        else:
            raise AttributeError('Tried to get unknown parameter {0}.'.format(param))
            # self.logger.warning('Tried to set unknown parameter %s. Blatantly ignored.', parameter)
            # should this raise an exception instead?

        unt = self.default_units.get(param, None)
        if unt is not None:
            if isinstance(target_list[param], np.matrix):
                unit_matrix = np.eye(len(target_list[param])) * _resolve_unit(unt, self._get_matching_unit)
                return target_list[param] * unit_matrix
            else:
                return target_list[param] * _resolve_unit(unt, self._get_matching_unit)
        else:
            return target_list[param]

    def _get_matching_unit(self, matchobj, indx=None):
        if indx is None:
            return self[matchobj.group(2)]
        else:
            return self[matchobj.group(2)][indx]

    # def _load_default_parameters_xml(self):
    #     """
    #     Loads the default parameter dictionary from the file specified in odsswr.conf. The default dictionary
    #     determines also what type the parameters should be.
    #     """
    #     self.editedParams = []  # reset
    #
    #     default_comp = Xmlc('default', self.toe)
    #
    #     self.params = default_comp.load()
    #     self.params_types = default_comp.typedic()

    def _load_default_parameters(self):
        """
        Loads the default parameter dictionary from the file specified in odsswr.conf. The default dictionary
        determines also what type the parameters should be.
        """
        self.editedParams = []  # reset

        for elname, el in default_comp.items():
            if el['type'] == self.toe:
                self.params = copy.deepcopy(el['properties'])
                self.default_units = copy.deepcopy(el['units'])
                try:
                    self._associated = copy.deepcopy(el['associated'])
                except KeyError:
                    self._associated = {}
                break
        else:
            raise TypeError('Could not  find a suitable reference object for {0} in the default file ("{1}")'
                            .format(self.toe, default_entities_path))

        self.params_types_raw = {k: type(v) for k, v in list(self.params.items()) + list(self._associated.items())}

    def _dump(self):
        """
        Prints all the parameters names and values. Debug use.
        """
        for parameter, value in self.params.items():
            print('<' + '{0} {1} = {2}'.format(type(value), parameter, value)[7:])

    def jsonize(self, all_params=False, flatten_mtx=True, using=default_entities_path):
        super_dikt = {'type': self.toe, 'name': self.name}
        if not self.isnamed():
            super_dikt['term_perm'] = self.term_perm

        pls_flat = {}
        pls_mtx = {}

        for parameter, value in self.params.items():
            if not all_params:
                if parameter not in self.editedParams:
                    continue
            if isinstance(value, np.matrix):
                pls_flat[parameter] = value.tolist()
            else:
                pls_flat[parameter] = value
            pls_mtx[parameter] = value

        if using is not None:
            pr_cmp = {'properties': pls_mtx, 'type': self.toe}
            compare = load_dictionary_json(using)
            for k, v in compare.items():
                try:
                    np.testing.assert_equal(pr_cmp, v)
                    dicky = {'type': self.toe,
                             'path': using,
                             'name': k}
                    if not self.isnamed():
                        dicky['term_perm'] = self.term_perm
                    return dicky
                except AssertionError:
                    pass

        if flatten_mtx:
            super_dikt['properties'] = pls_flat
        else:
            super_dikt['properties'] = pls_mtx

        super_dikt['depends'] = {}
        for depend_parameter, value in self._associated.items():
            super_dikt['depends'][depend_parameter] = value

        return super_dikt

    # this method has to be overridden by descendants who need something different from a pure parameter dump
    def fcs(self, **hookup):
        """
        Generates the string for creation of the object in OPENDSS. Needs to be passed the parameters proper not of the
        object itself, but of the object as being part of a circuit, such as, but not only, name and topological
        connections.
        """
        assert hookup == {}
        s1 = 'New ' + self.toe.split('_')[0] + '.' + self.name  # _ splitting to allow name personalization outside dss

        s2 = ''
        for parameter in self.editedParams:  # printing of non-default parameters only was preferred for better
            # readability of the returned string
            s2 = s2 + ' ' + parameter + '=' + _odssrep(self[parameter])
        return s1 + s2


# named entities are the abstract ones (linecodes...) rather than physical (lines, capacitors...). In other words,
# named entities are those that can exist and be instantiated outside and above the existence of a circuit.
# This is because there's conceptually no point in duplicating the same linecode, but instead it's conceptually correct
# to duplicate, where necessary, the same phisical items (storages...). So the names of the phisical items are, more
# appropriately, memorized into the graph, and their name is therefore tied to a circuit instance.

class _NamedDSSentity(_DSSentity):

    @classmethod
    def isnamed(cls):
        return True

    def __init__(self, name, xml_rep=None, **kwargs):
        super().__init__(xml_rep, **kwargs)
        self.name = name

    def fcs(self, **hookup):
        return super().fcs(**hookup)

# -------------------------------------------------------------
# LOADSHAPE CLASS WITH XML LOADER
# -------------------------------------------------------------


class CsvLoadshape:
    """
    Allows to specify a Loadshape that refers to a CSV file. Requires a path.
    The name of the loadshape will be the same as the file basename.
    Automatically recognizes if header is present or not.

    :param csv_path: the csv file path.
    :type csv_path: str
    :param column_scheme: A dictionary of one or more int:<'hour'|'mult'|'qmult'> couples that associate the column with one of the hour, mult and qmult properties.
    :type column_scheme: dict
    :param interval: the number of base units (see below) that constitute a step in the data. If unspecified or None the first column of the csv will be used as the vector of times, that can be, in general, non uniformly spaced.
    :type interval: float
    :param hour_fraction: the fraction of 1 HOUR that constitutes the base unit. Defaults to 60.0, that is 1/60th of an hour, that is, a minute.
    :type hour_fraction: float
    :param use_actual: setting this to False indicates that you want to use the data values in the csv to rescale the base value of the DSS object, rather than using the values directly.
    :type use_actual: bool
    :param npts: the number of points to load from the csv. If None, all the lines in the csv will be loaded in the loadshape. If greater than the number of lines in the csv, the lines will be tiled from the beginning. Please note that, since the automatic determination of the number of datapoints requires opening the csv and counting the rows, specifying this parameter, when possible, will grant a speed-up, especially for big files and/or multiple imports.
    :type npts: int

    """

    def __init__(self, csv_path: str, column_scheme: dict, interval=_pint_qty_type, use_actual=True, npts=None):

        if column_scheme == {}:
            raise ValueError('Empty column scheme')

        for kay in column_scheme.keys():
            if kay not in ('hour', 'mult', 'qmult'):
                raise ValueError('Unrecognized column key %s', kay)

        self._data = None
        self.column_scheme = column_scheme
        self.csv_path = os.path.abspath(csv_path)
        self.name = str(os.path.basename(csv_path)).split('.')[0]
        self.use_actual = use_actual
        if interval is None:
            self.intkey = 'interval'
            self.interval = 0.0
        else:  # using different properties just for cosmetic purposes
            if interval < 300 * um.s:
                self.intkey = 'sinterval'
                self.interval = interval.to(um.s).magnitude
            elif interval < 180 * um.min:
                self.intkey = 'minterval'
                self.interval = interval.to(um.min).magnitude
            else:
                self.intkey = 'interval'
                self.interval = interval.to(um.h).magnitude

        # auto-header recognition
        head = next(csv.reader(open(csv_path)))
        if all([_is_numeric_data(hitem) or _is_timestamp(hitem) for hitem in head]):
            self.header_string = 'No'
            self.shift = 0
        else:
            self.header_string = 'Yes'
            self.shift = 1

        # auto-row counting if no npts is passed
        if npts is None:  # if npts is not specified, automatic row counting is performed
            fo = csv.reader(open(csv_path))
            self.npts = str(sum([1 for row in fo]) - self.shift)  # -shift is for the header
        else:
            assert isinstance(npts, int)
            self.npts = str(npts)

        # auto-metadata recognizing
        row = next(csv.reader(open(csv_path)))  # the second row always has data in it
        ncol = len(row)

    @staticmethod
    def isnamed():
        return True

    @lru_cache()
    def get_at_row(self, row_no):

        if self._data is None:
            if self.shift == 0:
                pd_header = None
            elif self.shift == 1:
                pd_header = 0
            else:
                raise ValueError('Shift was neither 0 nor 1')
            self._data = read_csv(self.csv_path, header=pd_header)

        pcol = self.column_scheme.get('mult')
        qcol = self.column_scheme.get('qmult', None)

        mult = self._data.iloc[row_no, pcol-1]
        if qcol is not None:
            qmult = self._data.iloc[row_no, qcol-1]
        else:
            qmult = None

        return mult, qmult

    @lru_cache()
    def get_at_hour(self, hour):
        if not self.has_hour():
            raise AttributeError('Loadshape is defined by steps, not timestamps')
        hcol = self.column_scheme['hour']

        # hour should be sorted
        index = int(self._data[hcol].searchsorted(hour)) - 1

        return self.get_at_row(index)

    def has_hour(self):
        return 'hour' in self.column_scheme.keys()

    @property
    def uses_actual(self):
        return self.use_actual

    def fcs(self, **hookup):
        s = "New loadshape." + self.name + " npts=" + str(self.npts)

        s += self.intkey + "=" + str(self.interval)

        s += " Useactual=" + str(self.use_actual)

        for qty, ncol in self.column_scheme.items():
            s += " " + qty + "=(file=\"" + self.csv_path + "\", Column=" + str(ncol) + ", Header=" + self.header_string + ")"

        return s


class Loadshape:
    """Loadshape definition.
     DEPRECATED IN FAVOR OF CSVLOADSHAPE"""

    # diz = {'h': 'Interval',
    #        'm': 'mInterval',
    #        's': 'sInterval'}

    def __init__(self, name, p_shape_or_path=None, q_shape=None, interval_hours=1.0,
                 useactual=False):
        # yes, I know that this polimorphism is horrid, yet I did it

        assert isinstance(p_shape_or_path, (list, str))
        assert isinstance(name, str)
        self.name = name

        if isinstance(p_shape_or_path, list):  # if powers are specified, perform computed initialization
            p_shape = p_shape_or_path
            assert all([isinstance(item, float) for item in p_shape])

            self.p_mult = p_shape
            if q_shape is not None:  # if q_shape isn't specified, the same values of p.mult are used
                assert all([isinstance(item, float) for item in q_shape])
                self.q_mult = q_shape
            assert len(self.q_mult) == len(self.p_mult)

            if len([interval_hours]) != 1:  # a list was passed
                assert len(self.p_mult) == len([interval_hours])
                self.interval = 0.0  # in this case, interval has to be set to 0
                self.hour = interval_hours

            else:
                self.interval = interval_hours  # in this case, only interval is passed
                self.hour = None

            self.useactual = useactual

        else:  # if a path is specified, perform loading
            path = p_shape_or_path
            self.load(path, name)

    def load(self, path, name):
        # todo: for big loadshapes, it'd be useful to save the loadshape as a "stupid" csv and then writing mult=<file>

        try:
            parsedcontent = Xe.parse(path)
        except Xe.ParseError:
            raise
        except FileNotFoundError:
            raise

        root = parsedcontent.getroot()

        ls = root.find(".//loadshape[@name='{0}']".format(name))

        if ls is None:
            raise KeyError("Unable to find loadshape in the specified file")

        mult = ls.find('mult')
        interval = float(ls.find('interval').text)
        raw_useactual = ls.find('useactual').text
        if raw_useactual is None:
            self.useactual = False
        else:
            assert ls.find('useactual').text.lower() in ('true', 'false')
            if ls.find('useactual').text.lower() == 'true':
                self.useactual = True
            else:
                self.useactual = False

        # without checks, so an incomplete or empty curve raises exception
        self.p_mult = [float(pnt.find('p').text) for pnt in mult.iter('point')]

        raw_q_mult = [pnt.find('q').text for pnt in mult.iter('point')]
        if all(isinstance(a, str) for a in raw_q_mult):
            self.q_mult = list(map(lambda x: float(x), raw_q_mult))
        elif not (any(isinstance(a, str) for a in raw_q_mult)):
            self.q_mult = None
        else:
            raise ValueError("Data for the reactive power curve are incomplete")

        raw_hours = [pnt.find('t').text for pnt in mult.iter('point')]
        if all(isinstance(a, str) for a in raw_hours):
            self.hour = list(map(lambda x: float(x), raw_hours))
            self.interval = 0.0
        elif not (any(list(isinstance(a, str) for a in raw_hours))):
            self.hour = None
            assert isinstance(interval, float)
            assert interval > 0.0
            self.interval = interval
        else:
            raise ValueError("Data for the hours curve are incomplete")

    def npts(self):
        return len(self.p_mult)

    def fcs(self):
        s = 'New loadshape.' + self.name + ' '
        s += 'npts= ' + str(self.npts) + ' '
        s += 'interval= ' + str(self.interval) + ' '
        s += 'useactual= ' + str(self.useactual) + ' '

        shave = lambda array: '(' + str(array).replace(',', ' ').replace('[', '').replace(']', '') + ')'
        s += 'mult=' + shave(self.p_mult) + ' '

        if self.q_mult is not None:
            s += 'qmult=' + shave(self.q_mult) + ' '

        if self.hour is not None:
            s += 'hour=' + shave(self.hour)

        return s

    def plot(self):
        """Plots the loadshape"""
        # todo: add labels, correct ticks, title, etc
        # extend array head for plot correctness
        x = [0.0] + self.hour
        x = np.cumsum(x)
        y1 = [self.p_mult[0]] + self.p_mult
        if self.q_mult is None:
            y2 = y1
        else:
            y2 = [self.q_mult[0]] + self.q_mult

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.spines['bottom'].set_color('grey')
        ax.spines['top'].set_color('grey')
        ax.spines['right'].set_color('grey')
        ax.spines['left'].set_color('grey')
        ax.tick_params(axis='x', colors='#eeeeee')
        ax.tick_params(axis='y', colors='#eeeeee')
        ax.set_facecolor('#303030')
        ax.set_title('Loadshape [{0}]'.format(self.name), color='white')
        fig.patch.set_facecolor('#303030')
        fig.set_facecolor('#303030')
        plt.step(x, y1, label='Active load', color='#39ff14')
        plt.step(x, y2, label='Reactive load', color='#4cbb17', linestyle='--')
        plt.xlabel('Hours', color='white')
        plt.ylabel('Nominal power multiplier', color='white')
        plt.legend()
        plt.grid(True, 'both', linestyle='--', color='grey')
        plt.minorticks_on()
        plt.show()


# SUPPORT OBJECTS
# -------------------------------------------------------------

# note: these classes, along with the Line class, implement only an appropriate part of the underlying OpenDSS
# class properties, removing overlap between them.

class LineCode_S(_NamedDSSentity):
    """Contains a linecode defined with symmetrical components R0,R1,C0,C1,X0,X1.

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | nphases              | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | r1                   | float       | 0.05800                    |
    +----------------------+-------------+----------------------------+
    | x1                   | float       | 0.12060                    |
    +----------------------+-------------+----------------------------+
    | r0                   | float       | 0.17840                    |
    +----------------------+-------------+----------------------------+
    | x0                   | float       | 0.40470                    |
    +----------------------+-------------+----------------------------+
    | C1                   | float       | 3.40000                    |
    +----------------------+-------------+----------------------------+
    | C0                   | float       | 1.60000                    |
    +----------------------+-------------+----------------------------+
    | units                | string      | none                       |
    +----------------------+-------------+----------------------------+
    | baseFreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | normamps             | float       | 400.0                      |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 600.0                      |
    +----------------------+-------------+----------------------------+
    | faultrate            | float       | 0.1                        |
    +----------------------+-------------+----------------------------+
    | %perm                | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | repair               | float       | 3.0                        |
    +----------------------+-------------+----------------------------+
    | Kron                 | string      | N                          |
    +----------------------+-------------+----------------------------+
    | Rg                   | float       | 0.01805                    |
    +----------------------+-------------+----------------------------+
    | Xg                   | float       | 0.15508                    |
    +----------------------+-------------+----------------------------+
    | rho                  | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | neutral              | int         | 3                          |
    +----------------------+-------------+----------------------------+
    """
    pass


class LineCode_A(_NamedDSSentity):
    """Contains a linecode defined with generally asymmetrical components, rmatrix, cmatrix, xmatrix.

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | nphases              | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | units                | string      | none                       |
    +----------------------+-------------+----------------------------+
    | rmatrix              | floatmatrix | 0.09813333,0.04013333,.... |
    +----------------------+-------------+----------------------------+
    | xmatrix              | floatmatrix | 0.21530000,0.09470000,...  |
    +----------------------+-------------+----------------------------+
    | cmatrix              | floatmatrix | 2.80000000,-0.60000000,... |
    +----------------------+-------------+----------------------------+
    | baseFreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | normamps             | float       | 400.0                      |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 600.0                      |
    +----------------------+-------------+----------------------------+
    | faultrate            | float       | 0.1                        |
    +----------------------+-------------+----------------------------+
    | %perm                | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | repair               | float       | 3.0                        |
    +----------------------+-------------+----------------------------+
    | Kron                 | string      | N                          |
    +----------------------+-------------+----------------------------+
    | Rg                   | float       | 0.01805                    |
    +----------------------+-------------+----------------------------+
    | Xg                   | float       | 0.15508                    |
    +----------------------+-------------+----------------------------+
    | rho                  | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | neutral              | int         | 3                          |
    +----------------------+-------------+----------------------------+
    """
    pass


class WireData(_NamedDSSentity):
    """Contains overhead line conductor electro-physical constants and geometric parameters.

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | rdc                  | float       | -1                         |
    +----------------------+-------------+----------------------------+
    | rac                  | float       | -1                         |
    +----------------------+-------------+----------------------------+
    | runits               | string      | none                       |
    +----------------------+-------------+----------------------------+
    | gmrac                | float       | -1                         |
    +----------------------+-------------+----------------------------+
    | gmrunits             | string      | none                       |
    +----------------------+-------------+----------------------------+
    | radunits             | string      | none                       |
    +----------------------+-------------+----------------------------+
    | diam                 | float       | -2                         |
    +----------------------+-------------+----------------------------+
    """
    pass


class CNData(_NamedDSSentity):
    """Contains concentric neutral underground conductor electro-physical constant and geometric parameters

    +----------------------+-------------+----------------+
    | Property name        | type        | Default value  |
    +======================+=============+================+
    | k                    | int         | -1             |
    +----------------------+-------------+----------------+
    | DiaStrand            | float       | -1             |
    +----------------------+-------------+----------------+
    | GmrStrand            | float       |                |
    +----------------------+-------------+----------------+
    | Rstrand              | float       | -1             |
    +----------------------+-------------+----------------+
    | EpsR                 | float       |                |
    +----------------------+-------------+----------------+
    | InsLayer             | float       | -1             |
    +----------------------+-------------+----------------+
    | DiaIns               | float       |                |
    +----------------------+-------------+----------------+
    | DiaCable             | float       | -1             |
    +----------------------+-------------+----------------+
    | Rdc                  | float       | -1             |
    +----------------------+-------------+----------------+
    | Rac                  | float       | -2             |
    +----------------------+-------------+----------------+
    | Runits               | string      |                |
    +----------------------+-------------+----------------+
    | GMRac                | float       |                |
    +----------------------+-------------+----------------+
    | GMRunits             | string      |                |
    +----------------------+-------------+----------------+
    | radunits             | string      |                |
    +----------------------+-------------+----------------+
    | diam                 | float       |                |
    +----------------------+-------------+----------------+
    | k                    | float       | 2.3            |
    +----------------------+-------------+----------------+
    """
    pass


class TSData(_NamedDSSentity):
    """Contains tape shielded underground conductor electro-physical constant and geometric parameters

    +----------------------+-------------+---------------+
    | Property name        | type        | Default value |
    +======================+=============+===============+
    | DiaShield            | float       | -1            |
    +----------------------+-------------+---------------+
    | TapeLayer            | float       | -1            |
    +----------------------+-------------+---------------+
    | TapeLap              | float       |               |
    +----------------------+-------------+---------------+
    | EpsR                 | float       | -1            |
    +----------------------+-------------+---------------+
    | InsLayer             | float       | 0             |
    +----------------------+-------------+---------------+
    | DiaIns               | float       | -1            |
    +----------------------+-------------+---------------+
    | DiaCable             | float       |               |
    +----------------------+-------------+---------------+
    | Rdc                  | float       | -1            |
    +----------------------+-------------+---------------+
    | Rac                  | float       | -1            |
    +----------------------+-------------+---------------+
    | Runits               | string      | -2            |
    +----------------------+-------------+---------------+
    | GMRac                | float       |               |
    +----------------------+-------------+---------------+
    | GMRunits             | string      |               |
    +----------------------+-------------+---------------+
    | radunits             | string      |               |
    +----------------------+-------------+---------------+
    | normamps             | float       |               |
    +----------------------+-------------+---------------+
    | emergamps            | float       |               |
    +----------------------+-------------+---------------+
    | diam                 | float       |               |
    +----------------------+-------------+---------------+
    """
    pass


# Linegeom introduces a key difference in the storage of the data. Instead of listing the conductors in sequence
# (cond=1, prop1= value11, prop2=value21......propn1 = valuen1; cond=2, prop1 = value12, prop2 = value22.....)
# here the single conductor properties are stored in arrays whose index is ncond-1. The string forming function
# is then overridden in order to produce the declarations in opendss format.


class _LineGeometry(_NamedDSSentity):
    def __init__(self, name, xml_rep=None, **kwargs):
        super().__init__(name, xml_rep, **kwargs)
        ncond = self['nconds']
        self.specialparams = (self.wiretype, 'x', 'h', 'units')
        # for p in self.specialparams:
        #     if isinstance(self.params[p], np.matrix):
        #         assert self.params[p].size == ncond
        #     else:
        #         assert len(self.params[p]) == ncond

    def fcs(self, **hookup):
        s1 = 'New ' + self.toe.split('_')[0] + '.' + self.name

        s2 = ''
        for parameter in [p for p in self.editedParams if p not in self.specialparams]:
            s2 += ' ' + parameter + '=' + _odssrep(self[parameter])

        for ind in range(0, self['nconds']):
            s2 += '\n~ cond={0} '.format(ind + 1)
            for parameter in self.specialparams:
                if isinstance(self[parameter], _pint_qty_type):
                    true_param = self[parameter].magnitude
                else:
                    true_param = self[parameter]

                if isinstance(true_param, np.matrix):
                    idx = 0, ind  # matricial indicization necessary
                else:
                    idx = ind
                s2 += str(parameter) + '=' + str(true_param[idx]) + ' '

        return s1 + s2


class LineGeometry_O(_LineGeometry):
    """Line Geometry OVERHEAD

    +----------------------+-------------+----------------+
    | Property name        | type        | Default value  |
    +======================+=============+================+
    | nconds               | int         | 1              |
    +----------------------+-------------+----------------+
    | nphases              | int         | 1              |
    +----------------------+-------------+----------------+
    | wire                 | stringarray |                |
    +----------------------+-------------+----------------+
    | x                    | floatarray  | 0              |
    +----------------------+-------------+----------------+
    | h                    | floatarray  | 0              |
    +----------------------+-------------+----------------+
    | units                | stringarray | mt             |
    +----------------------+-------------+----------------+
    | reduce               | string      |                |
    +----------------------+-------------+----------------+
    | spacing              | string      |                |
    +----------------------+-------------+----------------+
    """

    def __init__(self, name, xml_rep=None, **kwargs):
        self.wiretype = 'wire'
        super().__init__(name, xml_rep, **kwargs)


class LineGeometry_T(_LineGeometry):
    """Line Geometry WITH TAPE SHIELDED CABLE

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | nconds               | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | nphases              | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | tscable              | stringarray |                            |
    +----------------------+-------------+----------------------------+
    | x                    | floatarray  | 0                          |
    +----------------------+-------------+----------------------------+
    | h                    | floatarray  | 0                          |
    +----------------------+-------------+----------------------------+
    | units                | stringarray | mt                         |
    +----------------------+-------------+----------------------------+
    | reduce               | string      |                            |
    +----------------------+-------------+----------------------------+
    | spacing              | string      |                            |
    +----------------------+-------------+----------------------------+
    """

    def __init__(self, name, xml_rep=None, **kwargs):
        self.wiretype = 'tscable'
        super().__init__(name, xml_rep, **kwargs)


class LineGeometry_C(_LineGeometry):
    """Line Geometry WITH CONCENTRIC NEUTRAL CABLE

    +----------------------+-------------+----------------+
    | Property name        | type        | Default value  |
    +======================+=============+================+
    | nconds               | int         | 1              |
    +----------------------+-------------+----------------+
    | nphases              | int         | 1              |
    +----------------------+-------------+----------------+
    | cncable              | stringarray |                |
    +----------------------+-------------+----------------+
    | x                    | floatarray  | 0              |
    +----------------------+-------------+----------------+
    | h                    | floatarray  | 0              |
    +----------------------+-------------+----------------+
    | units                | stringarray | mt             |
    +----------------------+-------------+----------------+
    | reduce               | string      |                |
    +----------------------+-------------+----------------+
    | spacing              | string      |                |
    +----------------------+-------------+----------------+
    """

    def __init__(self, name, xml_rep=None, **kwargs):
        self.wiretype = 'cncable'
        super().__init__(name, xml_rep, **kwargs)


class Curve:

    # todo implement csv and direct data polimorphism
    # todo port loadshape as curve

    _datadict = {'xycurve': ('xarray', 'yarray'),
                 'tshape': ('temp',)
    }

    def __init__(self, curve_type, name, data, interval=None, interval_unit='m'):
        self.name = name
        if curve_type not in ['xycurve', 'tshape']:
            raise ImportError('Unrecognized curve type')
        self.type = curve_type
        self.interval = interval
        self._dict = OrderedDict({'x': None, 'y': None, 'z': None})
        self.array_names = self._datadict[curve_type]

        assert interval_unit in ('h', 'm', 's')
        if interval_unit == 'h':
            string_interval_unit = ''
        else:
            string_interval_unit = interval_unit
        self.interval_string = string_interval_unit + 'interval'

        if isinstance(data, str) and data.endswith('.mat') or isinstance(data, dict):

            if isinstance(data, str):
                file_path = data
                mat = sio.loadmat(file_path, chars_as_strings=True)
            else:
                mat = data

            if curve_type in ['tshape']:
                assert self.interval is not None

            try:
                self._dict['x'] = mat['x']
                assert isinstance(self._dict['x'], np.ndarray)
            except KeyError:
                raise

            try:
                self._dict['y'] = mat['y']
                assert isinstance(self._dict['y'], np.ndarray)
                assert len(self._dict['x']) == len(self._dict['y'])
            except KeyError:
                self._dict['y'] = None

            try:
                self._dict['z'] = mat['z']
                assert isinstance(self._dict['z'], np.ndarray)
                assert len(self._dict['x']) == len(self._dict['z'])
            except KeyError:
                self._dict['z'] = None

            self.npts = len(self._dict['x'])

        elif isinstance(data, str) and data[-4:] == '.csv':

            csv_name = data
            self.npts = len(list(csv.reader(open(data)))) - 1

            for idx, nt in enumerate(self._dict.keys()):
                if idx >= len(self._datadict[curve_type]):
                    break
                self._dict[nt] = '(file=' + csv_name + ',Column=' + str(idx+1) + ',Header=Yes)'
        else:
            raise ImportError('Unrecognized data input. Please provide a .mat or .csv valid path or direct data.')

    @property
    def x(self):
        return self._dict['x']

    @property
    def y(self):
        return self._dict['y']

    @property
    def z(self):
        return self._dict['z']

    def isnamed(self):
        return True

    def fcs(self):
        s = 'New ' + self.type + '.' + self.name + ' ' + 'npts=' + str(self.npts) + ' '

        if self.interval is not None:
            s += self.interval_string + str(self.interval) + ' '

        for idx, array in enumerate([self.x, self.y, self.z]):
            if array is not None:
                s += self.array_names[idx] + '=' + str(array).replace('[[', '[').replace(']]', ']') + ' '

        return s


class PtCurve(Curve):
    pass


class EffCurve(Curve):
    pass

# -------------------------------------------------------------
# CIRCUIT ELEMENTS
# -------------------------------------------------------------
class _CircuitElement(_DSSentity):

    def __init__(self, xml_rep=None, **kwargs):
        super().__init__(xml_rep, **kwargs)

        # after loading from xml the line if appropriate, I go on with the kwargs in order to allow further on-the-fly
        # editing

        # self.current = None
        # self.voltage = None
        self.logger = None

        self._associated_elements = None

    def aka(self, name):
        try:
            assert self.name == ''
        except AssertionError:
            raise AssertionError(r'Cannot alias a component with name != ""')
        cpy = self()
        cpy.name = name
        return cpy


class _CircuitElementNBus(_CircuitElement):
    _nbusesdict = {'line': 2,
                  'reactor': 2,
                  'capacitor': 2,
                  'fault': 2,
                  'vsource': 2,
                  'isource': 1,
                  'generator': 1,
                  'storage': 1,
                  'load': 1,
                  'dload': 2,
                  'transformer': 2,
                  'switch': 2,
                  'pvsystem': 1,
                  'fourq': 1}
    # Format:
    # New elem.name bus1='a' bus2='b' ...busx='z' prop1=val1 prop2=val2....

    def __init__(self, xml_rep=None, **kwargs):
        super().__init__(xml_rep, **kwargs)
        self.nbuses = self._nbusesdict[self.toe]  # todo broken for transformers

    def fcs(self, **hookup):

        buses = hookup['buses']
        term_perm = hookup.get('terminals', None)

        s1 = 'New ' + self.toe.split('_')[0] + '.' + self.name  # _ splitting to allow name personalization outside dss
        s3 = ' '
        idox = 0
        for busno in range(1, self.nbuses + 1):
            if term_perm is not None:
                if isinstance(term_perm[(buses[busno - 1])], (tuple, int)):
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + _termrep(
                        term_perm[(buses[busno - 1])]) + ' '
                elif isinstance(term_perm[(buses[busno - 1])],
                                list):  # this happens when you specify more than one set of terminal connections at one bus
                    nthterminal = term_perm[(buses[busno - 1])][idox]
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + _termrep(nthterminal) + ' '
                    idox += 1
                elif term_perm[(buses[busno - 1])] is None:
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + ' '
            else:
                s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + ' '

        s2 = ''
        for parameter in self.editedParams:  # printing of non-default parameters only was preferred for better
            # readability of the returned string
            s2 = s2 + ' ' + parameter + '=' + _odssrep(self[parameter])
        return s1 + s3 + s2

    def _fes(self, params_to_indicate=None):
        # this function returns in the edit string all the requested params, even if not different from their default
        # values
        s1 = 'Edit ' + self.toe + '.' + self.name
        s2 = ''

        if params_to_indicate is None:
            params_selected = self.last_edited
        else:
            params_selected = params_to_indicate

        for parameter in params_selected:
            if not (parameter in self.params):
                self.logger.warning('Unknown parameter %s requested in edit string. Blatantly ignored.', parameter)
            else:
                s2 += ' ' + parameter + '=' + _odssrep(self[parameter])
        return s1 + s2


class Transformer(_DSSentity):  # remember that transformer is special, because has an arrayed mode of specifying
    # the memorized parameters
    """
    Transformer object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | windings             | int         | 2                          |
    +----------------------+-------------+----------------------------+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | conns                | stringarray | wye,wye                    |
    +----------------------+-------------+----------------------------+
    | kvs                  | floatarray  | 12.47,12.47                |
    +----------------------+-------------+----------------------------+
    | kvas                 | floatarray  | 1000.0,1000.0              |
    +----------------------+-------------+----------------------------+
    | taps                 | floatarray  | 1.000,1.000                |
    +----------------------+-------------+----------------------------+
    | %rs                  | floatarray  | 0.20,0.20                  |
    +----------------------+-------------+----------------------------+
    | xhl                  | float       | 7.000                      |
    +----------------------+-------------+----------------------------+
    | xht                  | float       | 35.000                     |
    +----------------------+-------------+----------------------------+
    | xlt                  | float       | 30.000                     |
    +----------------------+-------------+----------------------------+
    | x12                  | float       | 7.000                      |
    +----------------------+-------------+----------------------------+
    | x13                  | float       | 35.000                     |
    +----------------------+-------------+----------------------------+
    | x23                  | float       | 30.000                     |
    +----------------------+-------------+----------------------------+
    | xscmatrix            | floatmatrix | 7.00                       |
    +----------------------+-------------+----------------------------+
    | normmaxhkva          | float       | 1100.0                     |
    +----------------------+-------------+----------------------------+
    | emergmaxhkva         | float       | 1500.0                     |
    +----------------------+-------------+----------------------------+
    | thermal              | float       | 2.0                        |
    +----------------------+-------------+----------------------------+
    | n                    | float       | 0.8                        |
    +----------------------+-------------+----------------------------+
    | m                    | float       | 0.8                        |
    +----------------------+-------------+----------------------------+
    | flrise               | float       | 65.0                       |
    +----------------------+-------------+----------------------------+
    | hsrise               | float       | 15.0                       |
    +----------------------+-------------+----------------------------+
    | %loadloss            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | %noloadloss          | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | normhkva             | float       | 1100.0                     |
    +----------------------+-------------+----------------------------+
    | emerghkva            | float       | 1500.0                     |
    +----------------------+-------------+----------------------------+
    | sub                  | string      | n                          |
    +----------------------+-------------+----------------------------+
    | maxtap               | float       | 1.1                        |
    +----------------------+-------------+----------------------------+
    | mintap               | float       | 0.9                        |
    +----------------------+-------------+----------------------------+
    | numtaps              | int         | 32                         |
    +----------------------+-------------+----------------------------+
    | subname              | string      |                            |
    +----------------------+-------------+----------------------------+
    | %imag                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | ppm_antifloat        | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | bank                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | xfmrcode             | string      |                            |
    +----------------------+-------------+----------------------------+
    | xrconst              | string      | no                         |
    +----------------------+-------------+----------------------------+
    | x12                  | float       | 7.0                        |
    +----------------------+-------------+----------------------------+
    | x13                  | float       | 35.0                       |
    +----------------------+-------------+----------------------------+
    | x23                  | float       | 30.0                       |
    +----------------------+-------------+----------------------------+
    | leadlag              | string      | lag                        |
    +----------------------+-------------+----------------------------+
    | normamps             | float       | 50.929                     |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 69.449                     |
    +----------------------+-------------+----------------------------+
    | faultrate            | float       | 0.007                      |
    +----------------------+-------------+----------------------------+
    | %perm                | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | repair               | float       | 36.0                       |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+

    """

    def __init__(self, xml_rep=None, **kwargs):
        self.toe = 'transformer'
        super().__init__(xml_rep, **kwargs)

    def fcs(self, trname, windings, terminaldic):
        s2 = 'New transformer.' + trname

        for parameter in [x for x in self.editedParams if x not in ('conns', 'kvs', 'kvas', 'taps', '%rs')]:
            s2 += ' ' + parameter + '=' + _odssrep(self[parameter])

        s1 = ''
        for wdg, propdic in windings.items():
            if wdg[0] != trname + '_b':
                outbus = wdg[0]
            else:
                outbus = wdg[1]
            wgstr = '\n~ wdg=' + str(wdg[2] + 1) + ' '
            wgstr += 'bus=' + outbus + _termrep(terminaldic[wdg]) + ' '
            for propname, propvalue in propdic.items():
                wgstr += propname + '=' + _odssrep(propvalue) + ' '
            s1 += wgstr

        return s2 + s1


class Line(_CircuitElementNBus):
    """Line object. Data for the conductors must specified by a LineCode or a LineGeometry.

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | length               | float       | 1.000                      |
    +----------------------+-------------+----------------------------+
    | switch               | string      | false                      |
    +----------------------+-------------+----------------------------+
    | rg                   | float       | 0.01805                    |
    +----------------------+-------------+----------------------------+
    | xg                   | float       | 0.155081                   |
    +----------------------+-------------+----------------------------+
    | rho                  | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | units                | string      | none                       |
    +----------------------+-------------+----------------------------+
    | earthmodel           | string      | deri                       |
    +----------------------+-------------+----------------------------+
    | normamps             | float       | 400.0                      |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 600.0                      |
    +----------------------+-------------+----------------------------+
    """
    def __init__(self, xml_rep=None, **kwargs):
        # assert isinstance(data, (LineCode_A, LineCode_S, LineGeometry_C, LineGeometry_O, LineGeometry_T))
        super().__init__(xml_rep, **kwargs)

    def __call__(self, **kwargs):
        if 'data' in kwargs.keys():
            assert isinstance(kwargs['data'], (LineCode_S, LineCode_A, LineGeometry_T, LineGeometry_C, LineGeometry_O))
            self.data = kwargs['data']
            del kwargs['data']
        return super().__call__(**kwargs)

    # def fcs(self, **hookup):
    #     # todo use associated
    #     data_to_prop = {LineCode_A: 'linecode',
    #                     LineCode_S: 'linecode',
    #                     LineGeometry_T: 'geometry',
    #                     LineGeometry_O: 'geometry',
    #                     LineGeometry_C: 'geometry'}
    #
    #     #prop = data_to_prop[type(self.data)]
    #     #propname = self.data.name
    #
    #     base_string = super().fcs(**hookup)
    #     string = re.sub(' (?=bus1=)', ' ' + prop + '=' + propname + ' ', base_string, count=1)
    #
    #     return string

    # def jsonize(self, all_params=False, flatten_mtx=True, using=None):
    #     dk = super().jsonize(all_params, flatten_mtx, using)
    #     dk['depends'] = self['linecode'] + self['geometry']  # one of the two is blank
    #
    #     return dk

    def _setparameters(self, **kwargs):
        if 'data' in kwargs.keys():
            if isinstance(kwargs['data'], (LineCode_S, LineCode_A, LineGeometry_T, LineGeometry_C, LineGeometry_O)):
                self.data = kwargs['data']
                kwargs['data'] = kwargs['data'].name

        super()._setparameters(**kwargs)

    def __mul__(self, other):

        # todo pull

        if isinstance(other, (LineCode_A, LineCode_S, _LineGeometry)):
            super().__mul__(other)
            self['phases'] = other['nphases']
            return self
        else:
            raise KeyError



class Vsource(_CircuitElementNBus):
    """
    Voltage source

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | basekv               | float       | 115.0                      |
    +----------------------+-------------+----------------------------+
    | pu                   | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | angle                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | frequency            | float       | 60.0                       |
    +----------------------+-------------+----------------------------+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | mvasc3               | float       | 2000.0                     |
    +----------------------+-------------+----------------------------+
    | mvasc1               | float       | 2100.0                     |
    +----------------------+-------------+----------------------------+
    | x1r1                 | float       | 4.0                        |
    +----------------------+-------------+----------------------------+
    | x0r0                 | float       | 3.0                        |
    +----------------------+-------------+----------------------------+
    | isc3                 | float       | 10000.0                    |
    +----------------------+-------------+----------------------------+
    | isc1                 | float       | 10500.0                    |
    +----------------------+-------------+----------------------------+
    | r1                   | float       | 1.6038                     |
    +----------------------+-------------+----------------------------+
    | x1                   | float       | 6.4151                     |
    +----------------------+-------------+----------------------------+
    | r0                   | float       | 1.796                      |
    +----------------------+-------------+----------------------------+
    | x0                   | float       | 5.3881                     |
    +----------------------+-------------+----------------------------+
    | scantype             | string      | pos                        |
    +----------------------+-------------+----------------------------+
    | sequence             | string      | pos                        |
    +----------------------+-------------+----------------------------+
    | z1                   | floatarray  | 1.6037668,6.4150673        |
    +----------------------+-------------+----------------------------+
    | z0                   | floatarray  | 1.7960358,5.3881075        |
    +----------------------+-------------+----------------------------+
    | z2                   | floatarray  | 1.6037668,6.4150673        |
    +----------------------+-------------+----------------------------+
    | puz1                 | floatarray  | 0.012126781,0.048507125    |
    +----------------------+-------------+----------------------------+
    | puz0                 | floatarray  | 0.013580611,0.040741834    |
    +----------------------+-------------+----------------------------+
    | puz2                 | floatarray  | 0.012126781,0.048507125    |
    +----------------------+-------------+----------------------------+
    | basemva              | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | yearly               | string      |                            |
    +----------------------+-------------+----------------------------+
    | daily                | string      |                            |
    +----------------------+-------------+----------------------------+
    | duty                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | spectrum             | string      | defaultvsource             |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """
    pass


class Isource(_CircuitElementNBus):
    """
    Current source

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | amps                 | floatarray  | 0.0                        |
    +----------------------+-------------+----------------------------+
    | angle                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | frequency            | float       | 60.0                       |
    +----------------------+-------------+----------------------------+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | scantype             | string      | pos                        |
    +----------------------+-------------+----------------------------+
    | sequence             | string      | pos                        |
    +----------------------+-------------+----------------------------+
    | yearly               | string      |                            |
    +----------------------+-------------+----------------------------+
    | daily                | string      |                            |
    +----------------------+-------------+----------------------------+
    | duty                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | spectrum             | string      | default                    |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | enabled              | string      | true                       |
    +----------------------+-------------+----------------------------+
    """
    pass


class Fault(_CircuitElementNBus):
    """
    Fault object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | phases               | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | r                    | float       | 0.00                       |
    +----------------------+-------------+----------------------------+
    | %stddev              | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | ontime               | float       | 0.000                      |
    +----------------------+-------------+----------------------------+
    | temporary            | string      | no                         |
    +----------------------+-------------+----------------------------+
    | minamps              | float       | 5.0                        |
    +----------------------+-------------+----------------------------+
    | minamps              | float       | 5.0                        |
    +----------------------+-------------+----------------------------+
    | normamps             | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | faultrate            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | %perm                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | repair               | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """
    pass


class Capacitor(_CircuitElementNBus):
    """
    Capacitance

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | kvar                 | float       | 1200.0                     |
    +----------------------+-------------+----------------------------+
    | kv                   | float       | 12.470                     |
    +----------------------+-------------+----------------------------+
    | conn                 | string      | wye                        |
    +----------------------+-------------+----------------------------+
    | cuf                  | float       | 20.47                      |
    +----------------------+-------------+----------------------------+
    | r                    | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | xl                   | floatarray  | 0.0                        |
    +----------------------+-------------+----------------------------+
    | harm                 | intarray    | 0                          |
    +----------------------+-------------+----------------------------+
    | numsteps             | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | states               | intarray    | 1                          |
    +----------------------+-------------+----------------------------+
    | normamps             | float       | 75.0046059412345           |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 100.006141254979           |
    +----------------------+-------------+----------------------------+
    | faultrate            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | %perm                | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | repair               | float       | 3.0                        |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+

    """
    pass


class Reactor(_CircuitElementNBus):
    """Reactor object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | kvar                 | float       | 1200.0                     |
    +----------------------+-------------+----------------------------+
    | kv                   | float       | 12.47                      |
    +----------------------+-------------+----------------------------+
    | conn                 | string      | wye                        |
    +----------------------+-------------+----------------------------+
    | parallel             | string      | no                         |
    +----------------------+-------------+----------------------------+
    | r                    | float       | 0                          |
    +----------------------+-------------+----------------------------+
    | x                    | float       | 1555.009                   |
    +----------------------+-------------+----------------------------+
    | rp                   | float       | 0                          |
    +----------------------+-------------+----------------------------+
    | z1                   | floatarray  | 0,0                        |
    +----------------------+-------------+----------------------------+
    | z2                   | floatarray  | 0,0                        |
    +----------------------+-------------+----------------------------+
    | z0                   | floatarray  | 0,0                        |
    +----------------------+-------------+----------------------------+
    | z                    | floatarray  | 0.0,1555.009               |
    +----------------------+-------------+----------------------------+
    | rcurve               | string      |                            |
    +----------------------+-------------+----------------------------+
    | lcurve               | string      |                            |
    +----------------------+-------------+----------------------------+
    | lmh                  | float       | 4124.7895                  |
    +----------------------+-------------+----------------------------+
    | normamps             | float       | 5.0                        |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 6.0                        |
    +----------------------+-------------+----------------------------+
    | faultrate            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | %perm                | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | repair               | float       | 3.0                        |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
"""
    pass


class Generator(_CircuitElementNBus):
    """Generator object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | kv                   | float       | 12.47                      |
    +----------------------+-------------+----------------------------+
    | kw                   | float       | 1000.0                     |
    +----------------------+-------------+----------------------------+
    | pf                   | float       | 0.88                       |
    +----------------------+-------------+----------------------------+
    | kvar                 | float       | 60.0                       |
    +----------------------+-------------+----------------------------+
    | model                | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | vminpu               | float       | 0.90                       |
    +----------------------+-------------+----------------------------+
    | vmaxpu               | float       | 1.10                       |
    +----------------------+-------------+----------------------------+
    | yearly               | string      |                            |
    +----------------------+-------------+----------------------------+
    | daily                | string      |                            |
    +----------------------+-------------+----------------------------+
    | duty                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | dispmode             | string      | default                    |
    +----------------------+-------------+----------------------------+
    | dispvalue            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | conn                 | string      | wye                        |
    +----------------------+-------------+----------------------------+
    | rneut                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | xneut                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | status               | string      | variable                   |
    +----------------------+-------------+----------------------------+
    | class                | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | vpu                  | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | maxkvar              | float       | 120.0                      |
    +----------------------+-------------+----------------------------+
    | minkvar              | float       | -120.0                     |
    +----------------------+-------------+----------------------------+
    | pvfactor             | float       | 0.1                        |
    +----------------------+-------------+----------------------------+
    | forceon              | string      | no                         |
    +----------------------+-------------+----------------------------+
    | kva                  | float       | 1200.0                     |
    +----------------------+-------------+----------------------------+
    | mva                  | float       | 1.2                        |
    +----------------------+-------------+----------------------------+
    | xd                   | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | xdp                  | float       | 0.28                       |
    +----------------------+-------------+----------------------------+
    | xdpp                 | float       | 0.2                        |
    +----------------------+-------------+----------------------------+
    | h                    | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | d                    | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | usermodel            | string      |                            |
    +----------------------+-------------+----------------------------+
    | userdata             | string      |                            |
    +----------------------+-------------+----------------------------+
    | shaftmodel           | string      |                            |
    +----------------------+-------------+----------------------------+
    | shaftdata            | string      |                            |
    +----------------------+-------------+----------------------------+
    | dutystart            | float       | 0                          |
    +----------------------+-------------+----------------------------+
    | debugtrace           | string      | no                         |
    +----------------------+-------------+----------------------------+
    | balanced             | string      | no                         |
    +----------------------+-------------+----------------------------+
    | xrdp                 | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | spectrum             | string      | defaultgen                 |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """
    pass


class Storage(_CircuitElementNBus):
    """Storage object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | kv                   | float       | 12.47                      |
    +----------------------+-------------+----------------------------+
    | kw                   | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | pf                   | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | conn                 | string      | wye                        |
    +----------------------+-------------+----------------------------+
    | kvar                 | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | kva                  | float       | 25.0                       |
    +----------------------+-------------+----------------------------+
    | kwrated              | float       | 25.0                       |
    +----------------------+-------------+----------------------------+
    | kwhrated             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | kwhstored            | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | %stored              | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | %reserve             | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | state                | string      | idling                     |
    +----------------------+-------------+----------------------------+
    | %discharge           | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | %charge              | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | %effcharge           | float       | 90.0                       |
    +----------------------+-------------+----------------------------+
    | %effdischarge        | float       | 90.0                       |
    +----------------------+-------------+----------------------------+
    | %idlingkw            | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | %idlingkvar          | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | %r                   | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | %x                   | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | model                | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | vminpu               | float       | 0.9                        |
    +----------------------+-------------+----------------------------+
    | vmaxpu               | float       | 1.1                        |
    +----------------------+-------------+----------------------------+
    | balanced             | string      | no                         |
    +----------------------+-------------+----------------------------+
    | limitcurrent         | string      | no                         |
    +----------------------+-------------+----------------------------+
    | yearly               | string      |                            |
    +----------------------+-------------+----------------------------+
    | daily                | string      |                            |
    +----------------------+-------------+----------------------------+
    | duty                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | dispmode             | string      | default                    |
    +----------------------+-------------+----------------------------+
    | dischargetrigger     | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | chargetrigger        | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | timechargetrig       | float       | 2.0                        |
    +----------------------+-------------+----------------------------+
    | class                | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | dynadll              | string      |                            |
    +----------------------+-------------+----------------------------+
    | dynadata             | string      |                            |
    +----------------------+-------------+----------------------------+
    | usermodel            | string      |                            |
    +----------------------+-------------+----------------------------+
    | userdata             | string      |                            |
    +----------------------+-------------+----------------------------+
    | debugtrace           | string      | no                         |
    +----------------------+-------------+----------------------------+
    | spectrum             | string      |                            |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """
    pass


class Load(_CircuitElementNBus):
    """Load object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | kv                   | float       | 12.47                      |
    +----------------------+-------------+----------------------------+
    | kw                   | float       | 10.0                       |
    +----------------------+-------------+----------------------------+
    | pf                   | float       | 0.880                      |
    +----------------------+-------------+----------------------------+
    | model                | string      | 1                          |
    +----------------------+-------------+----------------------------+
    | yearly               | string      |                            |
    +----------------------+-------------+----------------------------+
    | daily                | string      |                            |
    +----------------------+-------------+----------------------------+
    | duty                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | growth               | float       |                            |
    +----------------------+-------------+----------------------------+
    | conn                 | string      | wye                        |
    +----------------------+-------------+----------------------------+
    | kvar                 | float       | 5.4                        |
    +----------------------+-------------+----------------------------+
    | rneut                | float       | -1.0                       |
    +----------------------+-------------+----------------------------+
    | xneut                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | status               | string      | variable                   |
    +----------------------+-------------+----------------------------+
    | class                | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | vminpu               | float       | 0.95                       |
    +----------------------+-------------+----------------------------+
    | vmaxpu               | float       | 1.05                       |
    +----------------------+-------------+----------------------------+
    | vminnorm             | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | vminemerg            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | xfkva                | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | allocationfactor     | float       | 0.500                      |
    +----------------------+-------------+----------------------------+
    | kva                  | float       | 11.4                       |
    +----------------------+-------------+----------------------------+
    | %mean                | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | %stddev              | float       | 10.0                       |
    +----------------------+-------------+----------------------------+
    | cvrwatts             | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | cvrvars              | float       | 2.0                        |
    +----------------------+-------------+----------------------------+
    | kwh                  | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | kwhdays              | float       | 30.0                       |
    +----------------------+-------------+----------------------------+
    | cfactor              | float       | 4.0                        |
    +----------------------+-------------+----------------------------+
    | cvrcurve             | string      |                            |
    +----------------------+-------------+----------------------------+
    | numcust              | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | zipv                 | floatarray  |                            |
    +----------------------+-------------+----------------------------+
    | %seriesrl            | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | relweight            | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | vlowpu               | float       | 0.5                        |
    +----------------------+-------------+----------------------------+
    | puxharm              | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | xrharm               | float       | 6.0                        |
    +----------------------+-------------+----------------------------+
    | spectrum             | string      | defaultload                |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """
    pass

    def __mul__(self, other):
        if isinstance(other, CsvLoadshape):
            self.edit(duty=other.name)
            return self
        else:
            raise ValueError('Object {0} is not addable to a Load'.format(str(other)))


class PvSystem(_CircuitElementNBus):

    """PvSystem object.

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | phases               | int         | 3                          |
    +----------------------+-------------+----------------------------+
    | kv                   | float       | 12.47                      |
    +----------------------+-------------+----------------------------+
    | irradiance           | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | Pmpp                 | float       | 500.0                      |
    +----------------------+-------------+----------------------------+
    | pctPmpp              | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | Temperature          | float       | 25.0                       |
    +----------------------+-------------+----------------------------+
    | pf                   | float       | 1.0                        |
    +----------------------+-------------+----------------------------+
    | conn                 | string      | wye                        |
    +----------------------+-------------+----------------------------+
    | kvar                 | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | kVA                  | float       | 500.0                      |
    +----------------------+-------------+----------------------------+
    | %Cutin               | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | %Cutout              | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | EffCurve             | float       |                            |
    +----------------------+-------------+----------------------------+
    | P-TCurve             | float       |                            |
    +----------------------+-------------+----------------------------+
    | %R                   | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | %X                   | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | model                | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | Vminpu               | float       | 0.9                        |
    +----------------------+-------------+----------------------------+
    | Vmaxpu               | float       | 1.1                        |
    +----------------------+-------------+----------------------------+
    | Balanced             | string      | No                         |
    +----------------------+-------------+----------------------------+
    | LimitCurrent         | string      | No                         |
    +----------------------+-------------+----------------------------+
    | yearly               | string      |                            |
    +----------------------+-------------+----------------------------+
    | daily                | string      |                            |
    +----------------------+-------------+----------------------------+
    | duty                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | Tyearly              | string      |                            |
    +----------------------+-------------+----------------------------+
    | Tdaily               | string      |                            |
    +----------------------+-------------+----------------------------+
    | Tduty                | string      |                            |
    +----------------------+-------------+----------------------------+
    | class                | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | UserModel            | string      |                            |
    +----------------------+-------------+----------------------------+
    | UserData             | string      |                            |
    +----------------------+-------------+----------------------------+
    | debugtrace           | string      | NO                         |
    +----------------------+-------------+----------------------------+
    | VarFollowInverter    | string      | No                         |
    +----------------------+-------------+----------------------------+
    | kvarLimit            | float       | 500.0                      |
    +----------------------+-------------+----------------------------+
    | DutyStart            | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | spectrum             | string      |                            |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """
    pass


# class Dload(_DSSentity):
#     """Distributed Load object. This is an original odsswr object, created from other base components of OPENDSS.
#
#     +----------------------+-------------+----------------------------+
#     | Property name        | type        | Default value              |
#     +======================+=============+============================+
#     | phases               | int         | 3                          |
#     +----------------------+-------------+----------------------------+
#     | mode                 | string      | lumped                     |
#     +----------------------+-------------+----------------------------+
#     | nsections            | int         | 12                         |
#     +----------------------+-------------+----------------------------+
#     """
#
#
#     def __init__(self, xml_rep=None, line=Line(), loads=None, **kwargs):
#         assert isinstance(line, Line)
#         if loads is not None:
#             for ld in loads.values():
#                 assert isinstance(ld, Load)
#         super().__init__(xml_rep, **kwargs)
#         self.line = line
#         self.loads = loads
#
#     def fcs(self, **hookup):
#
#         cname = self.name
#         buses = hookup['buses']
#         terminals = hookup.get('terminals', None)
#         # remember that loads is always a dictionary: {(term1,term2..) : Load}
#         totlength = self.line['length']
#
#         if self['mode'] == 'vdrop':
#             length_proportions = (1 / 2, 1 / 2)
#             load_proportions = (0, 1, 0)
#             modelization_type = 'VDROP EQUIVALENT'
#
#         elif self['mode'] == 'ploss':
#             length_proportions = (1 / 3, 2 / 3)
#             load_proportions = (0, 1, 0)
#             modelization_type = 'POWER LOSS EQUIVALENT'
#
#         elif self['mode'] == 'integral':
#             try:
#                 assert self['nsections'] is not None
#             except AssertionError:
#                 self.logger.error('distributed load %s was declared as integral but does not have nsections', cname)
#             ns = self['nsections']
#             length_proportions = [1 / ns] * ns
#             load_proportions = [0, *[1 / ns] * ns, 0]
#             modelization_type = 'INTEGRAL DISCRETIZED'
#
#         else:  # DEFAULTS TO LUMPED EXACT EQUIVALENT
#             if self['mode', 'none'] != 'lumped':  # so that if you explicit lumped mode, no info log is passed
#                 self.logger.info('distributed load %s defaulted to lumped exact equivalent mode')
#             length_proportions = (1 / 4, 3 / 4)
#             load_proportions = (0, 2 / 3, 1 / 3)
#             modelization_type = 'LUMPED EXACT EQUIVALENT'
#
#         n_segm = len(length_proportions)
#         midbuses = ['_DL_{0}#'.format(cname) + str(x) + '#' for x in range(1, n_segm)]
#         bus_sequence = [buses[0], *midbuses, buses[1]]
#
#         def get_terms(bus):
#             if terminals is not None:
#                 termd = {buses[1]: terminals[buses[1]]}
#                 return termd.get(bus, terminals[buses[0]])
#             else:
#                 return None
#
#         s = '! DISTRIBUTED LOAD (NAME:{0}) - {1}\n'.format(cname, modelization_type)
#
#         for ns in range(0, n_segm):
#             s += self.line(length=totlength * length_proportions[ns]).fcs([bus_sequence[ns], bus_sequence[ns + 1]],
#                                                                            '_L{0}_{1}'.format(ns + 1, cname),
#                                                                            {bus_sequence[ns]: get_terms(
#                                                                                bus_sequence[ns]),
#                                                                                bus_sequence[ns + 1]: get_terms(
#                                                                                    bus_sequence[ns + 1])}) + '\n'
#         for idxbus, bus in enumerate(bus_sequence):
#             if load_proportions[idxbus] != 0:
#                 for idxphase, (tm, ld) in enumerate(self.loads.items()):
#                     power = ld.params['kw']
#                     rpower = ld.params['kvar']
#                     s += ld(kw=power * load_proportions[idxbus], kvar=rpower * load_proportions[idxbus]).fcs(
#                         (bus,), '_V{0}({2})_{1}'.format(bus, cname, idxphase + 1), {bus: tm}) + '\n'
#
#         s += '! END DISTRIBUTED LOAD (NAME:{0})\n'.format(cname)
#
#         return s


class Switch(_CircuitElementNBus):
    """Switch object. The underlying model, if the switch is closed, behaves as a short, extremely conductive
    like trunk..

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | normamps             | float       | 400.0                      |
    +----------------------+-------------+----------------------------+
    | emergamps            | float       | 600.0                      |
    +----------------------+-------------+----------------------------+
    """
    def __init__(self, xml_rep=None, nphases=None, open=False):
        if open:
            self.disablestring = 'Disable = True'
        else:
            self.disablestring = ''
        super().__init__()
        self['phases'] = nphases

    def fcs(self, **hookup):
        buses = hookup['buses']
        s1 = '! SWITCH \nNew line.' + self.name  # _ splitting to allow name personalization outside dss
        s3 = ' '
        idox = 0
        for busno in range(1, self.nbuses + 1):
            if self.term_perm is not None:  # todo refactor with the new _termrep
                if isinstance(self.term_perm[(buses[busno - 1])], (tuple, int)):
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + _termrep(
                        self.term_perm[(buses[busno - 1])]) + ' '
                elif isinstance(self.term_perm[(buses[busno - 1])], list):
                    nthterminal = self.term_perm[(buses[busno - 1])][idox]
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + _termrep(nthterminal) + ' '
                    idox += 1
            else:
                s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + ' '

        s2 = ' switch=y length=0.1 r1=0.0001 r0=0.0001 x1=0.0 c1=0.0 x0=0.0 c0=0.0 ' + self.disablestring

        return s1 + s3 + s2

    def open(self):
        """Sets the switch OPEN."""
        self.disablestring = 'Disable = True'

    def close(self):
        """Sets the switch CLOSED."""
        self.disablestring = ''


# 4 quadrant, update capable component
class DecisionModel:
    @abstractmethod
    def decide_pq(self, oek, mynode):
        """Takes a graph and the node where the decision model has to interpret. Returns P(active power) and Q
        (reactive power), in that order, as a tuple. P and Q are not constrained in any way one to the other."""
        pass


class FourQ(Generator):

    @classmethod
    def isai(cls):
        return True

    def __init__(self, xml_rep=None, **kwargs):
        self._dm = None
        super().__init__(xml_rep, **kwargs)

    def update_pq(self, oek, mybus):
        assert self._dm is not None
        p, q = self._dm.decide_pq(oek, mybus)
        return p, q

    def fcs(self, **hookup):

        buses = hookup['buses']
        cname = self.name
        terminals = self.term_perm

        s1 = 'New ' + 'generator' + '.' + cname  # _ splitting to allow name personalization outside dss
        s3 = ' '
        idox = 0
        for busno in range(1, self.nbuses + 1):
            if terminals is not None:
                if isinstance(terminals[(buses[busno - 1])], (tuple, int)):
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + _termrep(
                        terminals[(buses[busno - 1])]) + ' '
                elif isinstance(terminals[(buses[busno - 1])],
                                list):  # this happens when you specify more than one set of terminal connections at one bus
                    nthterminal = terminals[(buses[busno - 1])][idox]
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + _termrep(nthterminal) + ' '
                    idox += 1
                elif terminals[(buses[busno - 1])] is None:
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + ' '
            else:
                s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + ' '

        s2 = ''
        for parameter in self.editedParams:  # printing of non-default parameters only was preferred for better
            # readability of the returned string
            s2 = s2 + ' ' + parameter + '=' + _odssrep(self[parameter])
        return s1 + s3 + s2

    def fus(self, oek, myname):
        mybus = oek[myname].topological['bus1']
        p, q = self.update_pq(oek, mybus)
        s = 'edit generator.' + myname + ' kw=' + str(p) + ' kvar=' + str(q)
        return s

    def define_dm(self, dm):
        assert isinstance(dm, DecisionModel)
        self._dm = dm

# ANCILLARY CLASSES
# -------------------------------------------------------------

class Capcontrol(_CircuitElement):
    pass

class Energymeter(_CircuitElement):
    pass

class Regcontrol(_CircuitElement):
    """Regcontrol object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | vreg                 | float       | 120.0                      |
    +----------------------+-------------+----------------------------+
    | band                 | float       | 3.0                        |
    +----------------------+-------------+----------------------------+
    | ptratio              | float       | 60.0                       |
    +----------------------+-------------+----------------------------+
    | ctprim               | float       | 300.0                      |
    +----------------------+-------------+----------------------------+
    | r                    | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | x                    | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | bus                  | float       |                            |
    +----------------------+-------------+----------------------------+
    | delay                | float       | 15.0                       |
    +----------------------+-------------+----------------------------+
    | reversible           | string      | no                         |
    +----------------------+-------------+----------------------------+
    | revvreg              | float       | 120.0                      |
    +----------------------+-------------+----------------------------+
    | revband              | float       | 3.0                        |
    +----------------------+-------------+----------------------------+
    | revr                 | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | revx                 | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | tapdelay             | float       | 2.0                        |
    +----------------------+-------------+----------------------------+
    | debugtrace           | string      | no                         |
    +----------------------+-------------+----------------------------+
    | maxtapchange         | int         | 16                         |
    +----------------------+-------------+----------------------------+
    | inversetime          | string      | no                         |
    +----------------------+-------------+----------------------------+
    | tapwinding           | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | vlimit               | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | ptphase              | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | revthreshold         | float       | 100.0                      |
    +----------------------+-------------+----------------------------+
    | revdelay             | float       | 60.0                       |
    +----------------------+-------------+----------------------------+
    | revneutral           | string      | no                         |
    +----------------------+-------------+----------------------------+
    | eventlog             | string      | yes                        |
    +----------------------+-------------+----------------------------+
    | remoteptratio        | float       | 60.0                       |
    +----------------------+-------------+----------------------------+
    | tapnum               | float       | 0                          |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """

    def __init__(self, xml_rep=None, **kwargs):
        super().__init__(xml_rep, **kwargs)

    def fcs(self, **hookup):
        transformer = hookup.get('trf')
        winding = hookup.get('wdg')
        s1 = ' transformer={0} winding={1}'.format(transformer, winding)

        s3 = super().fcs()

        return s3 + s1


class Monitor(_CircuitElement):
    """Monitor object. Several methods of odsswr.Circuit.py allow for exportation of the recordings of the monitors.

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | mode                 | int         | 0                          |
    +----------------------+-------------+----------------------------+
    | action               | string      |                            |
    +----------------------+-------------+----------------------------+
    | residual             | string      | no                         |
    +----------------------+-------------+----------------------------+
    | vipolar              | string      | yes                        |
    +----------------------+-------------+----------------------------+
    | ppolar               | string      | yes                        |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """

    def __init__(self, xml_rep=None, **kwargs):
        super().__init__(xml_rep, **kwargs)

    def fcs(self, **hookup):

        eltype = hookup.get('eltype')
        elname = hookup.get('elname')
        terminal = hookup.get('terminal')
        alias = hookup.get('alias')

        if alias is not None:
            name = alias
        else:
            name = 'mntr_' + eltype + '_' + elname

        return 'New monitor.' + name + ' element=' + eltype + '.' + elname + \
               ' terminal=' + str(terminal) + ' mode=' + str(self['mode']) + \
               ' vipolar=no ppolar=no'



class BusVoltageMonitor(_CircuitElement):
    def __init__(self, xml_rep=None, **kwargs):
        super().__init__(xml_rep, **kwargs)
        self['mode'] = 0

    def fcs(self, **hookup):

        busname = hookup['busname']

        l = Load(None, kw=0.0, kvar=0.0, model='1')
        l.name = 'bvm_' + busname

        loadstr = l.fcs(buses=(busname,))
        monstr = Monitor().fcs(eltype='load', elname='bvm_' + busname, terminal=1, alias='bvm_' + busname)

        return '!FICTITIOUS 0-KW LOAD FOR VOLTAGE MONITORING\n' + \
               loadstr + '\n' + \
               monstr + '\n' + \
               '!END FICTITIOUS LOAD'


class StorageController(_CircuitElement):
    """StorageController object

    +----------------------+-------------+----------------------------+
    | Property name        | type        | Default value              |
    +======================+=============+============================+
    | Element              | string      |                            |
    +----------------------+-------------+----------------------------+
    | Terminal             | int         | 1                          |
    +----------------------+-------------+----------------------------+
    | kWTarget             | float       | 8000.0                     |
    +----------------------+-------------+----------------------------+
    | %kWBand              | float       | 2.0                        |
    +----------------------+-------------+----------------------------+
    | PFTarget             | float       | 0.96                       |
    +----------------------+-------------+----------------------------+
    | PFBand               | float       | 0.04                       |
    +----------------------+-------------+----------------------------+
    | ElementList          | stringarray | s1                         |
    +----------------------+-------------+----------------------------+
    | Weights              | floatarray  | 1                          |
    +----------------------+-------------+----------------------------+
    | ModeDischarge        | string      | peakshave                  |
    +----------------------+-------------+----------------------------+
    | ModeCharge           | string      | Time                       |
    +----------------------+-------------+----------------------------+
    | TimeDischargeTrigger | int         | -1                         |
    +----------------------+-------------+----------------------------+
    | TimeChargeTrigger    | int         | 2                          |
    +----------------------+-------------+----------------------------+
    | %RatekW              | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | %Ratekvar            | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | %RateCharge          | float       | 20.0                       |
    +----------------------+-------------+----------------------------+
    | %Reserve             | float       | 25.0                       |
    +----------------------+-------------+----------------------------+
    | kWhTotal             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | kWTotal              | float       | 25.0                       |
    +----------------------+-------------+----------------------------+
    | kWhActual            | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    | kWActual             | float       | 5.0                        |
    +----------------------+-------------+----------------------------+
    | kWneed               | float       | 0.0                        |
    +----------------------+-------------+----------------------------+
    | %Participation       | float       |                            |
    +----------------------+-------------+----------------------------+
    | Yearly               | string      |                            |
    +----------------------+-------------+----------------------------+
    | Daily                | string      |                            |
    +----------------------+-------------+----------------------------+
    | Duty                 | string      |                            |
    +----------------------+-------------+----------------------------+
    | EventLog             | string      | No                         |
    +----------------------+-------------+----------------------------+
    | VarDispatch          | string      | No                         |
    +----------------------+-------------+----------------------------+
    | InhibitTime          | float       | 5.0                        |
    +----------------------+-------------+----------------------------+
    | Tup                  | float       | 0.25                       |
    +----------------------+-------------+----------------------------+
    | TFlat                | float       | 2.0                        |
    +----------------------+-------------+----------------------------+
    | Tdn                  | float       | 0.25                       |
    +----------------------+-------------+----------------------------+
    | kWThreshold          | float       | 6000.0                     |
    +----------------------+-------------+----------------------------+
    | basefreq             | float       | 50.0                       |
    +----------------------+-------------+----------------------------+
    """

    def __init__(self, xml_rep=None, **kwargs):
        super().__init__(xml_rep, **kwargs)

    def fcs(self, **hookup):

        eltype = hookup.get('eltype')
        elname = hookup.get('elname')
        terminal = hookup.get('terminal')
        alias = hookup.get('alias')

        if alias is not None:
            name = alias
        else:
            name = 'sc_' + eltype + '_' + elname

        s1 = 'New storagecontroller.' + name + ' element=' + eltype + '.' + elname + \
               ' terminal=' + str(terminal)

        s2 = ''
        for parameter in self.editedParams:  # printing of non-default parameters only was preferred for better
            # readability of the returned string
            s2 = s2 + ' ' + parameter + '=' + _odssrep(self[parameter])

        return s1 + s2


# ENGINE CLASS
# -------------------------------------------------------------
Element_c = namedtuple('Element', ['Name', 'Voltages', 'Currents', 'NumConductors', 'NumTerminals', 'NodeOrder'
                                   ])
Bus_c = namedtuple('Bus', ['Name', 'Voltages', 'VoltagesPU', 'Distance', 'NumNodes', 'kVBase'])
Monitor_c = namedtuple('Monitor', ['Name', 'Time', 'Header', 'NumChannels'])


class Element_x(Element_c):
    def __getattribute__(self, item):
        try:
            return getattr(super(), item)
        except AttributeError:
            try:
                currentname = odr.CktElement.Name()
                odr.Circuit.SetActiveElement(self.Name)
                fck = getattr(odr.CktElement, item)
                value = fck()
                odr.Circuit.SetActiveElement(currentname)
                return value
            except AttributeError:
                raise


class Engine:
    def __init__(self):
        assert odr.Basic.Start()
        self.dumps = {}

    @staticmethod
    def command(command_string):
        engine_command(command_string)

    @staticmethod
    def solve():
        engine_command('solve')
        pass

    def get_monitor(self, monitor_name):

        self._set_active_monitor(monitor_name)

        mn = Monitor_c(
            Name=monitor_name,
            Time=list(odr.Monitors.ByteStream()['t(sec)']/3600 + odr.Monitors.ByteStream()['hour']),
            Header=list(odr.Monitors.ByteStream().keys())[2:],
            NumChannels=len(list(odr.Monitors.ByteStream().keys())[2:])
        )

        return mn

    @staticmethod
    def _set_active_monitor(monitor_name):
        name_0 = odr.Monitors.Name()
        try:
            odr.Monitors.Name(monitor_name)
            assert odr.Monitors.Name().lower() == monitor_name.lower()
        except AssertionError:
            odr.Monitors.Name(name_0)
            assert odr.Monitors.Name().lower() == name_0.lower()
            raise KeyError('monitor name not found')

    def get_bus(self):
        pass

    def reset_time(self):
        pass

    def get_element(self):
        pass

    def squelch(self):
        pass

    def unsquelch(self):
        pass

    def init_mon_dump(self, path, n_points=None):

        if path in self.dumps.values():
            existing_index = [i for i, v in self.dumps.items() if v == path]
            assert len(existing_index) == 1
            return existing_index[0]
        else:
            try:
                index = max(self.dumps.keys()) + 1
            except ValueError:
                index = 0

        self.dumps[index] = path
        return index

    def dump_monitors(self, dumpid):
        for mn in odr.Monitors.AllNames():
            self.dump_monitor(mn, self.dumps[dumpid] + '\\' + mn + '.mat')

    def dump_monitor(self, monitor_name, path):

        self._set_active_monitor(monitor_name)
        bystr = odr.Monitors.ByteStream()
        raw_keyz_pairs = list(utils.aux_fcn.pairs(list(bystr.keys())[2:]))
        keyz = [re.sub(r'\W+', '', n1+n2) for n1, n2 in raw_keyz_pairs]
        out = DataFrame(columns=keyz)
        for k, (rk1, rk2) in zip(keyz, raw_keyz_pairs):
            out[k] = bystr[rk1] + 1j*bystr[rk2]
        # todo incremental dump
        sio.savemat(path, {'struct': out.to_dict('list')})

    @property
    def Buses(self):
        nm = odr.Circuit.AllBusNames()
        for busname in nm:
            odr.Circuit.SetActiveBus(busname)

            raw_volts = odr.Bus.Voltages()
            assert len(raw_volts) % 2 == 0
            volts = []
            for rp, ip in zip(raw_volts[0::2], raw_volts[1::2]):
                volts.append(rp + 1j*ip)

            raw_puvolts = odr.Bus.PuVoltage()
            assert len(raw_puvolts) % 2 == 0
            puvolts = []
            for rp, ip in zip(raw_puvolts[0::2], raw_puvolts[1::2]):
                puvolts.append(rp + 1j * ip)

            b = Bus_c(
                Name=odr.Bus.Name(),
                Voltages=volts,
                VoltagesPU=puvolts,
                Distance=odr.Bus.Distance(),
                NumNodes=odr.Bus.NumNodes(),
                kVBase=odr.Bus.kVBase())

            yield b

    @property
    def Lines(self):
        for el in self.Elements:
            if re.match('^Line\.', el.Name):
                yield el
            else:
                continue

    @property
    def Elements(self):
        nm = odr.Circuit.AllElementNames()
        for elname in nm:
            odr.Circuit.SetActiveElement(elname)

            raw_volts = odr.CktElement.Voltages()
            assert len(raw_volts) % 2 == 0
            volts = []
            for rp, ip in zip(raw_volts[0::2], raw_volts[1::2]):
                volts.append(rp + 1j*ip)

            raw_currs = odr.CktElement.Currents()
            assert len(raw_currs) % 2 == 0
            currs = []
            for rp, ip in zip(raw_currs[0::2], raw_currs[1::2]):
                currs.append(rp + 1j*ip)

            nc = odr.CktElement.NumConductors()
            nt = odr.CktElement.NumTerminals()

            el = Element_x(
                Name=elname,  # todo treat
                Voltages=np.reshape(volts, (nt, nc)),
                Currents=np.reshape(currs, (nt, nc)),
                NumConductors=nc,
                NumTerminals=nt,
                NodeOrder=odr.CktElement.NodeOrder())

            yield el

# MAIN FUNCTION FOR DEMONSTRATION AND TESTING
# -------------------------------------------------------------

def main():
    pass

if __name__ == "__main__":
    main()
