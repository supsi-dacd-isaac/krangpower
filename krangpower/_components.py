# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import copy
import csv
import hashlib
import json
import os.path
import re
import sys
import textwrap
from abc import abstractmethod
from collections import OrderedDict
from functools import lru_cache

import numpy as np
import scipy.io as sio
from pandas import read_csv

from ._aux_fcn import termrep, is_numeric_data, is_timestamp, matrix_from_json
from ._config_loader import PINT_QTY_TYPE, DEFAULT_ENTITIES_PATH, ASSOCIATION_TYPES_PATH, \
    UM, DEFAULT_COMP, DSSHELP, GLOBAL_PRECISION, TMP_PATH, MANDATORY_UNITS
from ._exceptions import AssociationError, TypeUnrecoverableError, RecoveryTargetError, TypeRecoveryError
from ._logging_init import mlog
from ._nxtable import NxTable

__all__ = ['CsvLoadshape', 'LineGeometry_C', 'LineGeometry_T', 'LineGeometry_O',
           'LineCode_A', 'LineCode_S', 'Line', 'WireData', 'CNData', 'TSData', 'Curve', 'PtCurve', 'EffCurve',
           'Vsource', 'dejsonize', 'SnpMatrix', 'load_entities', 'LineCode', 'LineGeometry',
           'Isource', 'Load', 'Transformer', 'Capacitor', 'Capcontrol', 'Regcontrol', 'Reactor',
           'Monitor', 'StorageController', 'Storage', 'PvSystem', 'Generator', 'FusAble']

_muldict = NxTable()
_fmtdict = NxTable()
_muldict.from_csv(ASSOCIATION_TYPES_PATH, ccol=3)
_fmtdict.from_csv(ASSOCIATION_TYPES_PATH, ccol=4)


def get_classmap():

    comp_module = sys.modules[__name__]
    classmap = {}
    for item in dir(comp_module):
        classmap[item.lower()] = getattr(comp_module, item)

    return classmap


# <editor-fold desc="AUX FUNCTIONS">

def _cpx(odss_tuple, nterm, ncond):
    """
    This function transforms the raw data for electric parameters (voltage, current...) in a suitable complex array

    :param odss_tuple: tuple of nphases*2 floats (returned by odsswr as couples of real, imag components, for each phase
        of each terminal)
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

    if isinstance(data_raw, PINT_QTY_TYPE):
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

    elif isinstance(data, SnpMatrix):

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


# </editor-fold>


# <editor-fold desc="AUX CLASSES">
class SnpMatrix(np.matrix):  # extends np.matrix, allowing to instantiate a symmetrical mtx by passing a tril string.
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


def dejsonize(obj_repr: dict):

    classmap = get_classmap()

    def propgetter(matchobj, indx=None):

        if indx is None:
            try:
                return obj_repr['properties'][matchobj.group(2)]
            except KeyError:
                return DEFAULT_COMP['default_' + obj_repr['type']]['properties'][matchobj.group(2)]
        else:
            try:
                return obj_repr['properties'][matchobj.group(2)][indx]
            except KeyError:
                return DEFAULT_COMP['default_' + obj_repr['type']]['properties'][matchobj.group(2)]['indx']

    # determines class
    elcls = classmap[obj_repr['type']]

    if 'path' in obj_repr.keys():
        with open(obj_repr['path'], 'r') as file:
            dik = json.load(file)
            obj_repr['properties'] = dik[obj_repr['name']]['properties']

    # restore matrices
    for prop, value in obj_repr['properties'].items():
        if isinstance(value, list):
            obj_repr['properties'][prop] = matrix_from_json(value)

    # add unit measure
    for prop in obj_repr['properties'].keys():

        if prop not in DEFAULT_COMP['default_' + obj_repr['type']]['units'].keys():
            continue
        else:
            try:
                value = obj_repr['units'][prop]
            except KeyError:
                mlog.debug('When djsonizing {0}.{1}, property {0} did not have a unit, so the default unit was assumed'.format(
                    obj_repr['type'], obj_repr['name'], prop))
                value = DEFAULT_COMP['default_' + obj_repr['type']]['units'][prop]

        if isinstance(obj_repr['properties'][prop], np.matrix):
            unit_matrix = np.eye(len(obj_repr['properties'][prop])) * _resolve_unit(value, propgetter)
            obj_repr['properties'][prop] = np.multiply(obj_repr['properties'][prop], unit_matrix)
        elif obj_repr['properties'][prop] is None:
            pass
        else:
            obj_repr['properties'][prop] *= _resolve_unit(value, propgetter)

    # add in the adjointed parameters
    obj_repr['properties'].update(obj_repr['depends'])

    # returns object
    if elcls.isnamed():
        return elcls(obj_repr['name'], **obj_repr['properties'])
    else:
        return elcls(**obj_repr['properties'])


def _type_recovery(value, target_type):
    try:
        if isinstance(value, int):
            assert target_type == float
            recovered_value = float(value)

        elif isinstance(value, (list, np.ndarray)):
            assert target_type == np.matrix
            recovered_value = np.matrix(value)
        elif isinstance(value, str):
            recovered_value = target_type(value)
        elif isinstance(value, np.complex):
            assert target_type == np.matrix
            recovered_value = np.matrix([np.real(value), np.imag(value)])
        else:
            raise TypeUnrecoverableError(type(value))
    except AssertionError:
        raise RecoveryTargetError(type(value), target_type)
    except ValueError:
        raise RecoveryTargetError(type(value), target_type)

    return recovered_value


def _resolve_unit(ustring: str, match_unit_getter):

    # resolution of self-referencing units, of the type "ohm / this.units" where "units" is a string parameter
    # of the _DSSObj
    if not isinstance(ustring, (list, tuple)):
        istring = [ustring]
    else:
        istring = ustring

    units = []

    # todo after unit dereferencing, this is no more necessary

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
        return [UM.parse_units(s) for s in units][0]
    else:
        return [UM.parse_units(s) for s in units]


def _get_help(config, cmp):
    section = cmp.eltype
    helpitems = config.items(section.upper().split('_')[0])
    help_str = ''
    basev = 90

    for item in helpitems:
        dname = item[0].split(') ')
        pname = dname[-1].lower().replace(' ', '')
        if pname in DEFAULT_COMP['default_' + section]['properties'].keys():
            mod = ''
            def_val = str(DEFAULT_COMP['default_' + section]['properties'][pname])
            if def_val == '':
                def_val = '<none>'
            def_value = '\n\t ' + '[default value=' + def_val + ']'
            try:
                cur_value = '\n\t ' + '[current value=' + str(cmp[pname].magnitude) + ']'
            except AttributeError:
                cur_value = '\n\t ' + '[current value=' + str(cmp[pname]) + ']'
        else:
            mod = '\u0336'
            def_value = ''
            cur_value = ''
        # hstr = mod.join(item[0].upper() + def_value + '\n\t' + item[1])

        wrapper = textwrap.TextWrapper(width=basev*(len(mod)+1), replace_whitespace=True)
        help_str += mod.join(item[0].upper() + def_value + cur_value + '\n\t ' +'\n\t '.join(wrapper.wrap(item[1])) + '\n')
    return help_str


# -------------------------------------------------------------
# GENERIC DSSENTITY AND ITS NAMED VERSION
# -------------------------------------------------------------
class FcsAble:
    def __init__(self):
        self._multiplied_objs = []

    @classmethod
    def isnamed(cls):
        """Returns True if the class is a named entity, like a LineCode, False otherwise."""
        return False

    @classmethod
    def isabove(cls):
        """Returns True if the class is an above-graph component, like a monitor or regulator, False otherwise."""
        return False

    @classmethod
    def isai(cls):
        """Returns true if the class is ai-enabled, False otherwise."""
        return False

    @abstractmethod
    def fcs(self, **hookup):
        pass


class FusAble(FcsAble):

    @classmethod
    def isai(cls):
        """Returns true if the class is ai-enabled, False otherwise."""
        return True

    @abstractmethod
    def fus(self, oek, myname):
        pass


class _DSSentity(FcsAble):
    def __init__(self, **parameters):
        super().__init__()
        self.term_perm = None
        self.name = ''
        self.toe = self.__class__.__name__.lower()
        self._params = None
        self._editedParams = None
        self._default_units = None

        self._load_default_parameters()

        self._setparameters(**parameters)

    def __str__(self):
        return 'krangpower.' + self.__class__.__name__ + '(' + self.name + ')'

    def __repr__(self):
        return '<' + self.__str__() + '(@' + str(hex(id(self))) + ')>'

    @property
    def eltype(self):
        return self.toe.split('_')[0]

    @property
    def fullname(self):
        """Name, in the form eltype.elname"""
        return (self.eltype + '.' + self.name).lower()

    def paramhelp(self):
        """Prints a cheatsheet for the object's parameters."""
        print('\nPARAMETERS HELP FOR {0} (get/set them with {0}[<param>])\n'.format(self.eltype))
        print(_get_help(DSSHELP, self))
        return None

    def aka(self, name):
        """Aliases the object."""
        # try:
        #     assert self.name == ''
        # except AssertionError:
        #     raise AssertionError(r'Cannot alias a component with name != ""')
        cpy = self()
        cpy.name = name
        return cpy

    def sf_deepcopy(self):
        """Returns a deep copy of the object. This method is to be used in place of the standard copy.deepcopy, that
        will raise an exception."""
        # the problems are caused by what is possibly inside _multiplied_objs, so we work around the issue by copying
        # the object without it and then reassigning the content.
        memo = self._multiplied_objs
        self._multiplied_objs = []
        cpy = copy.deepcopy(self)
        cpy._multiplied_objs = memo
        self._multiplied_objs = memo

        return cpy

    def __imul__(self, other):

        try:
            i2 = other.eltype
        except AttributeError:
            i2 = other[0].eltype  # happens for lists of components

        try:
            prop_to_set = _muldict[self.toe, i2]
            prop_fmt = _fmtdict[self.toe, i2]
        except KeyError:
            try:
                oname = other.name
            except AttributeError:
                oname = ''
            raise AssociationError(i2, oname, self.eltype, self.name)

        # support for setting ATTRIBUTES of the object, beginning with '.'
        if prop_to_set.startswith('.'):
            attr_to_set = prop_to_set[1:]
            if prop_fmt == 'self':  # this means that the value to be set is the object itself
                setattr(self, attr_to_set, other)
            else:
                setattr(self, attr_to_set, getattr(other, prop_fmt))

        # support for setting PROPERTIES of the object
        else:
            if isinstance(other, list):
                self[prop_to_set] = [getattr(o, prop_fmt) for o in other]
                self._multiplied_objs.extend(other)
            else:
                self[prop_to_set] = getattr(other, prop_fmt)
                self._multiplied_objs.append(other)

            # if we set a property, we also memorize object in this list, so, when declaring the object, we can
            # also automatically declare objects that were multiplied. this is not needed for attributes!
            # self._multiplied_objs.extend(other)

    def __mul__(self, other):
        """Associates the multiplier with this object according to the built-in association rules."""
        
        belf = self.sf_deepcopy()
        
        i1 = belf.toe
        try:
            i2 = other.eltype
        except AttributeError:
            i2 = other[0].eltype  # happens for lists of components

        try:
            prop_to_set = _muldict[i1, i2]
            prop_fmt = _fmtdict[i1, i2]
        except KeyError:
            try:
                oname = other.name
            except AttributeError:
                oname = ''
            raise AssociationError(i2, oname, belf.eltype, belf.name)

        # support for setting ATTRIBUTES of the object, beginning with '.'
        if prop_to_set.startswith('.'):
            attr_to_set = prop_to_set[1:]
            if prop_fmt == 'self':  # this means that the value to be set is the object itbelf
                setattr(belf, attr_to_set, other)
            else:
                setattr(belf, attr_to_set, getattr(other, prop_fmt))

        # support for setting PROPERTIES of the object
        else:
            if isinstance(other, list):
                belf[prop_to_set] = [getattr(o, prop_fmt) for o in other]
                belf._multiplied_objs.extend(other)
            else:
                belf[prop_to_set] = getattr(other, prop_fmt)
                belf._multiplied_objs.append(other)

            # if we set a property, we also memorize object in this list, so, when declaring the object, we can
            # also automatically declare objects that were multiplied. this is not needed for attributes!
            # belf._multiplied_objs.extend(other)

        return belf

    def __call__(self, **kwargs):
        # a call returns a copy of the object edited on the fly with the kwargs passed
        edited = copy.deepcopy(self)
        edited._setparameters(**kwargs)
        return edited

    # def __getattr__(self, item):
    #     # fallback getattr tries to call getparameters
    #     try:
    #         tp = self._params[item]
    #     except KeyError:
    #         raise AttributeError
    #
    #     return tp

    def __getitem__(self, item):
        """Gets a parameter for the object."""
        return self._getparameter(item)

    def __setitem__(self, key, value):
        """Sets a parameter for the object."""
        self._setparameters(**{key: value})

    def _setparameters(self, **kwargs):

        def unitsfirst(dicky):
            for k, v in dicky.items():
                if k in ('units', 'runits', 'gmrunits'):
                    # after unit dereferencing, this should only fire when loading from dss and in general when
                    # passing raw input
                    mlog.debug('The unitsfirst parameter setting routine detected an unit setting with {0}={1}.'
                               .format(k, str(v)))
                    yield k, v
                else:
                    continue

            for k, v in dicky.items():
                if k not in ('units', 'runits', 'gmrunits'):
                    yield k, v
                else:
                    continue

        for parameter_raw, value_raw in unitsfirst(kwargs):

            # check against direct setting of secured parameters
            if parameter_raw.lower() in self._ignored_params:
                raise ValueError('You cannot directly set property {0} for an object of type {1}. This property has '
                                 'to be set by associating an entity of the correct type.'
                                 .format(parameter_raw, self.toe))

            # patch for the 'pct in keyword' problem
            if re.match('pct', parameter_raw.lower()):
                parameter = parameter_raw.lower().replace('pct', '%')
            else:
                parameter = parameter_raw

            # select what dict regards the parameter: _params, _associated?
            if parameter.lower() in self._params.keys():
                target_list = self._params
                default_dict = 'properties'
            elif parameter.lower() in self._associated.keys():
                default_dict = 'associated'
                target_list = self._associated

            else:
                raise AttributeError('Tried to set unknown parameter {0}.'.format(parameter))

            # pint quantity check and conversion
            if hasattr(value_raw, 'magnitude'):  # isinstance pint does not work
                unt = _resolve_unit(self._default_units[parameter.lower()], self._get_prop_from_matchobj)
                if unt == UM.none:
                    raise ValueError('There is no unit for {0}. This should not happen, contact the dev.'
                                     .format(self.toe))
                value = value_raw.to(unt).magnitude
            else:
                value = value_raw

            # type-correctness check and, if failed, recovery attempt, of the raw value
            try:
                assert isinstance(value, (self.params_types_raw[parameter.lower()], type(None)))
            except AssertionError:
                target_type = self.params_types_raw[parameter.lower()]
                try:
                    value = _type_recovery(value, target_type)

                # if the value couldn't be salvaged, raise
                except TypeRecoveryError:
                    raise RecoveryTargetError(type(value),
                                              target_type,
                                              msg='parameter "{0}" is of type {1} instead of {2}'
                                                  ' and could not be converted'.format(parameter.lower(),
                                                                                type(value),
                                                                                target_type))

            if isinstance(value, np.matrix):
                test = np.array_equal(value, DEFAULT_COMP['default_' + self.toe][default_dict][parameter])
            elif isinstance(value, list):
                test = value == DEFAULT_COMP['default_' + self.toe][default_dict][parameter.lower()]
            else:
                test = value == DEFAULT_COMP['default_' + self.toe][default_dict][parameter.lower()]

            if test:
                mlog.debug('[{2}-{3}]Ignored setting {0} = {1} because identical to default'
                           .format(parameter, str(value), self.toe, self.name))
                continue

            # finally setting the parameter
            target_list[parameter.lower()] = value
            self._editedParams.append(parameter.lower())

    def _getparameter(self, param):

        if param.lower() in self._params.keys():
            target_list = self._params
        elif param.lower() in self._associated.keys():
            target_list = self._associated
        else:
            raise AttributeError('Tried to get unknown parameter {0}.'.format(param))
            # self.logger.warning('Tried to set unknown parameter %s. Blatantly ignored.', parameter)
            # should this raise an exception instead?

        unt = self._default_units.get(param, None)
        if unt is not None:
            if isinstance(target_list[param], np.matrix):
                unit_matrix = np.eye(len(target_list[param])) * _resolve_unit(unt, self._get_prop_from_matchobj)
                return target_list[param] * unit_matrix
            else:
                return target_list[param] * _resolve_unit(unt, self._get_prop_from_matchobj)
        else:
            return target_list[param]

    def _get_prop_from_matchobj(self, matchobj, indx=None):
        if indx is None:
            return self[matchobj.group(2)]
        else:
            return self[matchobj.group(2)][indx]

    def _load_default_parameters(self):
        """
        Loads the default parameter dictionary from the file specified in odsswr.conf. The default dictionary
        determines also what type the parameters should be.
        """
        self._editedParams = []  # reset

        for elname, el in DEFAULT_COMP.items():
            if el['type'] == self.toe:
                self._params = copy.deepcopy(el['properties'])
                self._default_units = copy.deepcopy(el['units'])
                try:
                    self._associated = copy.deepcopy(el['associated'])
                except KeyError:
                    self._associated = {}

                try:
                    self._ignored_params = copy.deepcopy(el['ignored'])
                except KeyError:
                    self._ignored_params = []

                break
        else:
            raise TypeError('Could not  find a suitable reference object for {0} in the default file ("{1}")'
                            .format(self.toe, DEFAULT_ENTITIES_PATH))

        self.params_types_raw = {k: type(v) for k, v in list(self._params.items()) + list(self._associated.items())}

    def jsonize(self, all_params=False, flatten_mtx=True):
        """Returns a dict that describes the element's parameters. It's compatible with the json I/O functions."""
        super_dikt = {'type': self.toe, 'name': self.name, 'units': {}}
        if not self.isnamed():
            super_dikt['term_perm'] = self.term_perm

        pls_flat = {}
        pls_mtx = {}

        for parameter, value in self._params.items():
            if not all_params:
                if parameter not in self._editedParams:
                    continue
            if isinstance(value, np.matrix):
                pls_flat[parameter] = np.round(value, GLOBAL_PRECISION).tolist()
            elif isinstance(value, float):
                pls_flat[parameter] = np.round(value, GLOBAL_PRECISION)
            else:
                pls_flat[parameter] = value
            pls_mtx[parameter] = value

        # all units in _params are in standard units. Therefore, we pass standard units.
        super_dikt['units'] = {k: v for k, v in DEFAULT_COMP['default_' + self.toe]['units'].items()
                               if k in self._editedParams}

        # if using is not None:
        #     pr_cmp = {'properties': pls_mtx, 'type': self.toe}
        #     compare = load_dictionary_json(using)
        #     for k, v in compare.items():
        #         try:
        #             np.testing.assert_equal(pr_cmp, v)
        #             dicky = {'type': self.toe,
        #                      'path': using,
        #                      'name': k}
        #             if not self.isnamed():
        #                 dicky['term_perm'] = self.term_perm
        #             return dicky
        #         except AssertionError:
        #             pass

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
        """Generates the opendss "new..." command for the object. Needs to be passed an appropriate set of keyword
         arguments to fill whatever topologic parameters are needed. It's meant primarily for internal use by
         krangpower itself."""
        assert hookup == {}
        s1 = 'New ' + self.toe.split('_')[0] + '.' + self.name  # _ splitting to allow name personalization outside dss

        s2 = ''

        if self.toe in MANDATORY_UNITS.keys():
            for muname, muvalue in MANDATORY_UNITS[self.toe].items():
                s2 = s2 + ' ' + muname + '=' + muvalue

        # we now declare the different-from-default parameters only for better readability of the returned string

        # we first declare the associated parameters, because (for opendss buggish behavior) they are required
        # at the beginning even if declared with their keyword
        for parameter in [p for p in self._editedParams if p in self._associated.keys()]:
            s2 = s2 + ' ' + parameter + '=' + _odssrep(self[parameter])

        for parameter in [p for p in self._editedParams if p not in self._associated.keys()]:
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

    def __init__(self, name='', **parameters):
        # can be instantiated with a blank name, but you won't be able to add it to a krang without aka-ing it!
        super().__init__(**parameters)
        self.name = name


# -------------------------------------------------------------
# LOADSHAPE CLASS
# -------------------------------------------------------------
class Loadshape(_DSSentity):
    pass


class CsvLoadshape(FcsAble):
    """
    Allows to specify a Loadshape that refers to a CSV file. Requires a path.
    The name of the loadshape will be the same as the file basename.
    Automatically recognizes if header is present or not.

    :param name: the loadshape name.
    :type name: str
    :param csv_path: the csv file path.
    :type csv_path: str
    :param column_scheme: A dictionary of one or more int:<'hour'|'mult'|'qmult'> couples that associate the column with
     one of the hour, mult and qmult properties.
    :type column_scheme: dict
    :param interval: length of step in the data as a pint unit with dimension [time]. If unspecified or None the first
     column of the csv will be used as the vector of times, that can be, in general, non uniformly spaced.
    :type interval: _pint_qty_type
    :param use_actual: setting this to False indicates that you want to use the data values in the csv to rescale the
     base value of the DSS object, rather than using the values directly.
    :type use_actual: bool
    :param npts: the number of points to load from the csv. If None, all the lines in the csv will be loaded in the
     loadshape. If greater than the number of lines in the csv, the lines will be tiled from the beginning. Please note
     that, since the automatic determination of the number of datapoints requires opening the csv and counting the rows,
     specifying this parameter, when possible, will grant a speed-up, especially for big files and/or multiple imports.
    :type npts: int
    """

    eltype = 'CsvLoadshape'

    def __init__(self, name='', csv_path=None, column_scheme=None, interval=None, use_actual=True, npts=None):
        super().__init__()
        if column_scheme == {}:
            raise ValueError('Empty column scheme')

        if column_scheme is None:
            column_scheme = {'mult': 1}

        for kay in column_scheme.keys():
            if kay not in ('hour', 'mult', 'qmult'):
                raise ValueError('Unrecognized column key %s', kay)

        self._data = None
        self.column_scheme = column_scheme

        # if csv_path is a fully qualified path, we load it; if it's just the name, we look for the file
        # inside TMP_PATH (where csvs are unpacked, for instance, by open_ckt); if nothing works, we raise
        if os.path.isfile(os.path.abspath(csv_path)):
            self.csv_path = os.path.abspath(csv_path)
        elif os.path.isfile(os.path.join(TMP_PATH, name + '.csv')):
            self.csv_path = os.path.join(TMP_PATH, name + '.csv')
        else:
            raise IOError('Could not find file {0}'.format(csv_path))

        if name == '':
            self.name = str(os.path.basename(self.csv_path)).split('.')[0]
        else:
            self.name = name
        self.use_actual = use_actual
        if isinstance(interval, PINT_QTY_TYPE):
            self.true_interval = interval
        else:
            self.true_interval = interval * UM.min

        if interval is None:
            self.intkey = 'interval'
            self.interval = 0.0
        else:  # using different properties just for cosmetic purposes
            if self.true_interval < (300 * UM.s):
                self.intkey = 'sinterval'
                self.interval = self.true_interval.to(UM.s).magnitude
            elif self.true_interval < (180 * UM.min):
                self.intkey = 'minterval'
                self.interval = self.true_interval.to(UM.min).magnitude
            else:
                self.intkey = 'interval'
                self.interval = self.true_interval.to(UM.h).magnitude

        # auto-header recognition
        head = next(csv.reader(open(self.csv_path)))
        if all([is_numeric_data(hitem) or is_timestamp(hitem) for hitem in head]):
            self.header_string = 'No'
            self.shift = 0
        else:
            self.header_string = 'Yes'
            self.shift = 1

        # auto-row counting if no npts is passed
        if npts is None:  # if npts is not specified, automatic row counting is performed
            fo = csv.reader(open(csv_path))
            self.npts = str(sum([1 for _ in fo]) - self.shift)  # -shift is for the header
        else:
            assert isinstance(npts, int)
            self.npts = str(npts)

        self._calchash()

        # # auto-metadata recognizing
        # row = next(csv.reader(open(csv_path)))  # the second row always has data in it
        # ncol = len(row)

    @staticmethod
    def isnamed():
        return True

    @staticmethod
    def isabove():
        return False

    @property
    def fullname(self):
        return 'csvloadshape.' + self.name

    def _calchash(self):
        # hash
        with open(self.csv_path, 'r') as csvfile:
            cnt = csvfile.read()
        self.hash = hashlib.md5(cnt.encode('utf-8')).hexdigest()
        return self.hash

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
        s = "New loadshape." + self.name + " npts=" + str(self.npts) + " "

        s += self.intkey + "=" + str(self.interval)

        s += " Useactual=" + str(self.use_actual)

        for qty, ncol in self.column_scheme.items():
            s += " " + qty + "=(file=\"" + self.csv_path + "\", Column=" + str(ncol) + ", Header=" + self.header_string + ")"

        return s

    def jsonize(self):

        super_dikt = {'type': 'csvloadshape', 'name': self.name, 'depends': {}, 'properties': {}, 'units': {}}
        super_dikt['properties']['csv_path'] = os.path.basename(self.csv_path)
        super_dikt['properties']['column_scheme'] = self.column_scheme
        super_dikt['properties']['npts'] = int(self.npts)
        super_dikt['properties']['use_actual'] = self.use_actual
        super_dikt['properties']['interval'] = self.true_interval.to('min').magnitude
        super_dikt['units']['interval'] = str(self.true_interval.units)
        super_dikt['hash'] = self._calchash()

        return super_dikt

    def peek(self):
        import matplotlib
        matplotlib.use('QT5agg')
        import matplotlib.pyplot as plt

        data = np.genfromtxt(self.csv_path, delimiter=',', skip_header=self.shift, names=['kw', 'kvar'])
        ics = np.asarray(list(range(int(self.npts)))) * self.true_interval.to('hour').magnitude
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(ics, data['kw'], color='k', label='kw')
        ax1.plot(ics, data['kvar'], color='r', label='kvar')
        ax1.set_title(self.name)
        plt.show()


# SUPPORT OBJECTS
# -------------------------------------------------------------

# note: these optisim, along with the Line class, implement only an appropriate part of the underlying OpenDSS
# class properties, removing overlap between them.

class LineCode(_NamedDSSentity):

    def jsonize(self, all_params=False, flatten_mtx=True):
        sudict = super().jsonize(all_params, flatten_mtx)
        if not np.any(np.isnan(sudict['properties']['r0'])):
            for prop in ('rmatrix', 'xmatrix', 'cmatrix'):
                try:
                    del sudict['properties'][prop]
                except KeyError:
                    continue
        else:
            for prop in ('r0', 'x0', 'c0', 'r1', 'x1', 'c1'):
                try:
                    del sudict['properties'][prop]
                except KeyError:
                    continue

        return sudict


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


class LineGeometry(_NamedDSSentity):
    def __init__(self, name, **parameters):
        self.wiretype = None
        super().__init__(name, **parameters)

        # self.specialparams
        # for p in self.specialparams:
        #     if isinstance(self._params[p], np.matrix):
        #         assert self._params[p].size == ncond
        #     else:
        #         assert len(self._params[p]) == ncond

    @property
    def specialparams(self):
        return self.wiretype, 'x', 'h'

    def fcs(self, **hookup):

        if self.wiretype is None:
            raise ValueError('No wire was specified for {}'.format(self.fullname))

        s1 = 'New ' + self.toe.split('_')[0] + '.' + self.name

        s2 = ''
        for parameter in [p for p in self._editedParams if p not in self.specialparams]:
            s2 += ' ' + parameter + '=' + _odssrep(self[parameter])

        for ind in range(0, self['nconds']):
            s2 += '\n~ cond={0} units=m '.format(ind + 1)
            for parameter in self.specialparams:
                if isinstance(self[parameter], PINT_QTY_TYPE):
                    true_param = self[parameter].magnitude
                else:
                    true_param = self[parameter]

                if isinstance(true_param, np.matrix):
                    idx = 0, ind  # matricial indicization necessary
                else:
                    idx = ind
                # try:
                #     # for fullnames of the wires
                #     s2 += str(parameter) + '=' + str(true_param[idx]).split('.')[1] + ' '
                # except IndexError:
                s2 += str(parameter) + '=' + str(true_param[idx]) + ' '
        return s1 + s2

    def _setparameters(self, **kwargs):
        pts = set(kwargs.keys())
        ccs = {'wire', 'cncable', 'tscable'}
        wire_set = pts.intersection(ccs)
        if len(wire_set) == 0:
            pass
        elif len(wire_set) == 1:
            if self.wiretype is not None:
                raise ValueError('The cable type for linegeometry.{} is already set!'.format(self.name))
            self.wiretype = list(wire_set)[0]
        else:
            raise ValueError('Tried to set more than one cabletype')

        super()._setparameters(**kwargs)


class LineGeometry_O(LineGeometry):
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

    def __init__(self, name, **parameters):
        super().__init__(name, **parameters)
        self.wiretype = 'wire'


class LineGeometry_T(LineGeometry):
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

    def __init__(self, name, **parameters):
        super().__init__(name, **parameters)
        self.wiretype = 'tscable'


class LineGeometry_C(LineGeometry):
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

    def __init__(self, name, **parameters):
        super().__init__(name, **parameters)
        self.wiretype = 'cncable'


class Curve(FcsAble):

    # todo implement csv and direct data polimorphism
    # todo port loadshape as curve

    _datadict = {'xycurve': ('xarray', 'yarray'),
                 'tshape': ('temp',)
    }

    def __init__(self, name, curve_type, data, interval=None):
        super().__init__()
        self.name = name
        if curve_type not in ['xycurve', 'tshape']:
            raise ImportError('Unrecognized curve type')
        self.type = curve_type
        self.interval = interval
        self._dict = OrderedDict({'x': None, 'y': None, 'z': None})
        self.array_names = self._datadict[curve_type]

        if isinstance(data, str) and data.endswith('.mat') or isinstance(data, dict):

            if isinstance(data, str):
                file_path = data
                mat = sio.loadmat(file_path, chars_as_strings=True)
            else:
                mat = data

            self._mat = mat

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

        elif isinstance(data, np.matrix):  # a matrix in which each row is x,y,z
            self._dict['x'] = np.asarray(data[0, :])

            try:
                self._dict['y'] = np.asarray(data[1, :])
            except IndexError:
                self._dict['y'] = None

            try:
                self._dict['z'] = np.asarray(data[2, :])
            except IndexError:
                self._dict['z'] = None

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

    @staticmethod
    def isnamed():
        return True

    def jsonize(self, all_params=False, flatten_mtx=True):

        super_dikt = {'type': 'curve', 'name': self.name, 'depends': {}, 'properties': {}, 'units': {}}
        super_dikt['properties']['curve_type'] = self.type
        # super_dikt['properties']['name'] = self.name
        try:
            super_dikt['properties']['interval'] = self.interval.magnitude
            super_dikt['units']['interval'] = str(self.interval.units)
        except AttributeError:  # happens for None!
            super_dikt['properties']['interval'] = self.interval

        if flatten_mtx:
            super_dikt['properties']['data'] = [self._dict[x].tolist()[0] for x in self._dict.keys() if self._dict[x] is not None]
        else:
            super_dikt['properties']['data'] = [self._dict[x] for x in self._dict.keys() if self._dict[x] is not None]

        return super_dikt

    def fcs(self):
        s = 'New ' + self.type + '.' + self.name + ' '

        if self.interval is not None:
            s += 'sinterval=' + str(self.interval.to('s').magnitude) + ' '

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

    def __init__(self, **parameters):
        super().__init__(**parameters)
        # after loading from xml the line if appropriate, I go on with the kwargs in order to allow further on-the-fly
        # editing


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
                   'switch': 2,
                   'pvsystem': 1,
                   'fourq': 1}
    # Format:
    # New elem.name bus1='a' bus2='b' ...busx='z' prop1=val1 prop2=val2....

    def __init__(self, **parameters):
        super().__init__(**parameters)
        self.nbuses = self._nbusesdict[self.toe]

    def fcs(self, **hookup):

        buses = hookup['buses']
        term_perm = hookup.get('terminals', None)

        s1 = 'New ' + self.toe.split('_')[0] + '.' + self.name  # _ splitting to allow name personalization outside dss

        s4 = ''
        if self.toe in MANDATORY_UNITS.keys():
            for muname, muvalue in MANDATORY_UNITS[self.toe].items():
                s4 = s4 + ' ' + muname + '=' + muvalue

        s3 = ' '
        idox = 0
        for busno in range(1, self.nbuses + 1):
            if term_perm is not None:
                if isinstance(term_perm[busno - 1], (tuple, int)):
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + termrep(
                        term_perm[busno - 1]) + ' '
                elif isinstance(term_perm[busno - 1], list):
                    # this happens when you specify more than one set of terminal connections at one bus
                    nthterminal = term_perm[busno - 1][idox]
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + termrep(nthterminal) + ' '
                    idox += 1
                elif term_perm[busno - 1] is None:
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + ' '
            else:
                s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + ' '

        s2 = ''
        for parameter in self._editedParams:  # printing of non-default parameters only was preferred for better
            # readability of the returned string
            s2 = s2 + ' ' + parameter + '=' + _odssrep(self[parameter])
        return s1 + s4 + s3 + s2


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

    def __init__(self, **parameters):
        self.toe = 'transformer'
        super().__init__(**parameters)
        self.specialparams = ('conns', 'kvs', 'kvas', 'taps', '%rs')

    def fcs(self, **hookup):

        buses = hookup['buses']
        termdic = hookup.get('terminals', None)

        s1 = 'New ' + self.toe.split('_')[0] + '.' + self.name

        s2 = ''
        for parameter in [p for p in self._editedParams if p not in self.specialparams]:
            s2 += ' ' + parameter + '=' + _odssrep(self[parameter])

        for ind in range(0, self['windings']):
            s2 += '\n~ wdg={0} bus={1}{2}'.format(ind + 1, buses[ind], termrep(termdic[ind])) + ' '
            for parameter in self.specialparams:
                if isinstance(self[parameter], PINT_QTY_TYPE):
                    true_param = self[parameter].magnitude
                else:
                    true_param = self[parameter]

                if isinstance(true_param, np.matrix):
                    idx = 0, ind  # matricial indicization necessary
                else:
                    idx = ind
                s2 += str(parameter).strip('s') + '=' + str(true_param[idx]) + ' '

        return s1 + s2

    def aka(self, name):
        try:
            assert self.name == ''
        except AssertionError:
            raise ValueError(r'Cannot alias a component with name != ""')
        cpy = self()
        cpy.name = name
        return cpy


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
    def __init__(self, **parameters):
        # assert isinstance(data, (LineCode_A, LineCode_S, LineGeometry_C, LineGeometry_O, LineGeometry_T))
        super().__init__(**parameters)

    def __call__(self, **kwargs):
        if 'data' in kwargs.keys():
            assert isinstance(kwargs['data'], (LineCode_S, LineCode_A, LineGeometry_T, LineGeometry_C, LineGeometry_O))
            self.data = kwargs['data']
            del kwargs['data']
        return super().__call__(**kwargs)

    def _setparameters(self, **kwargs):
        if 'data' in kwargs.keys():
            if isinstance(kwargs['data'], (LineCode_S, LineCode_A, LineGeometry_T, LineGeometry_C, LineGeometry_O)):
                self.data = kwargs['data']
                kwargs['data'] = kwargs['data'].name

        super()._setparameters(**kwargs)

    def __mul__(self, other):

        bsmul = super().__mul__(other)  # first, so it throws errors if necessary
        self['phases'] = other['nphases']
        return bsmul


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
    def __init__(self, nphases=None, is_open=False):
        if is_open:
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
                if isinstance(self.term_perm[busno - 1], (tuple, int)):
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + termrep(
                        self.term_perm[busno - 1]) + ' '
                elif isinstance(self.term_perm[busno - 1], list):
                    nthterminal = self.term_perm[busno - 1][idox]
                    s3 += 'bus' + str(busno) + '=' + str(buses[busno - 1]) + termrep(nthterminal) + ' '
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


# ANCILLARY CLASSES
# -------------------------------------------------------------
class _AboveCircuitElement(_DSSentity):
    @classmethod
    def isabove(cls):
        return True

    def __init__(self, name='', **parameters):
        # can be instantiated with a blank name, but you won't be able to add it to a krang without aka-ing it!
        super().__init__(**parameters)
        self.name = name


class Capcontrol(_AboveCircuitElement):
    pass


class Energymeter(_AboveCircuitElement):
    pass


class Regcontrol(_AboveCircuitElement):
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


class Monitor(_AboveCircuitElement):
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
    pass


class StorageController(_AboveCircuitElement):
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
    pass


def load_entities(path):

    if isinstance(path, str):
        with open(path, 'r') as file:
            dik = json.load(file)

    elif isinstance(path, dict):
        dik = path

    else:
        raise TypeError

    dicky = {}

    for entity_name in dik:

        eltype, elname = entity_name.split('.')
        if '_' in eltype:
            dep_name = ''.join(eltype.split('_')[:-1]) + '.' + elname
        else:
            dep_name = entity_name
        dicky[dep_name] = dejsonize(dik[entity_name])

    return dicky


# MAIN FUNCTION FOR DEMONSTRATION AND TESTING
# -------------------------------------------------------------
def main():
    pass


if __name__ == "__main__":
    main()

