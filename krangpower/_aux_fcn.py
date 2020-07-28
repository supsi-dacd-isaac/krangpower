# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import difflib
import io
import re
import textwrap
from sys import modules
from functools import lru_cache
import json
import hashlib
import canonicaljson
import numpy as np
from dateutil.parser import parse as dateparse


def fingerprint_file(path):
    with open(path, 'r') as file:
        md = json.load(file)

    return hashlib.md5(canonicaljson.encode_canonical_json(md)).hexdigest()


def get_help_out(config, section):
    helpitems = config.items(section.upper().split('_')[0])
    help_str = ''
    basev = 90
    for item in helpitems:
        hstr = item[0].upper() + ':  ' + item[1]
        help_str += '\n'+'\n\t '.join(textwrap.wrap(hstr, basev))
    return help_str


def diff_dicts(original: dict, new: dict, context_lines=1):

    oj = io.StringIO()
    json.dump(original, oj, sort_keys=True, indent=2)
    original_json = oj.getvalue()

    nj = io.StringIO()
    json.dump(new, nj, sort_keys=True, indent=2)
    new_json = nj.getvalue()

    diff = difflib.context_diff(original_json.splitlines(1),
                                new_json.splitlines(1),
                                n=context_lines, fromfile='original', tofile='new')
    err = ''.join(diff)

    return err


def ebus(bus: str, nt: int):
    return bus + '_' + str(nt)


def pairs(iterable):
    itr = iter(iterable)
    while True:
        try:
            yield next(itr), next(itr)
        except StopIteration:
            raise


def lower(item):

    if hasattr(item, 'lower'):
            return item.lower()
    elif hasattr(item, '__iter__'):
        try:
            return [s.lower() for s in item]
        except AttributeError:
            raise AttributeError('Not all the items contained in the argument have a "lower" method')
    else:
        raise AttributeError('The argument does not have a "lower" method')


def pairwise(iterable):
    # "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)


def matrix_from_json(value):

    def desym(lol):
        size = len(lol)
        dsm = np.zeros([size, size])
        for r in range(size):
            for c in range(r+1):
                dsm[r, c] = lol[r][c]
                dsm[c, r] = lol[r][c]

        return dsm

    if isinstance(value[0], str):
        return value
    else:
        try_mtx = np.asarray(value)
        if try_mtx.dtype == 'object':
            return desym(value)
        else:
            return try_mtx


@lru_cache(8)
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
                dik[entity]['properties'][prop] = matrix_from_json(value)

    return dik


def bus_resolve(bus_descriptor: str):
    """
    >>> bus_resolve('bus2.3.1.2')
    ('bus2', (3, 1, 2))
    >>> bus_resolve('bus2.33.14.12323.2.3.3')
    ('bus2', (33, 14, 12323, 2, 3, 3))
    """

    bus_descriptor.replace('bus.', '')
    tkns = bus_descriptor.split('.')

    bus = tkns[0]
    terminals = tuple(int(x) for x in tkns[1:])

    return bus, terminals


def is_timestamp(item):
    try:
        dateparse(item)
    except ValueError:
        return False
    return True


def termrep(terminals):
    """
    This function takes a terminal collection (represented by a tuple of ints) and returns a representation that can be
    cat to a bus name in order to form a full bus-terminal qualification according to the odsswr syntax.

    >>> termrep(1,3,2)
    '.1.3.2'

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
                s += str(t) + '.'
            return s[0:-1]  # shaves final dot
        except TypeError:
            return '.' + str(terminals)


def is_numeric_data(item):
    return re.fullmatch('([0-9]|\,|\.| )*', item) is not None