import json
import textwrap
from sys import modules
from functools import lru_cache

import numpy as np


def get_help_out(config, section):
    helpitems = config.items(section.upper().split('_')[0])
    help_str = ''
    basev = 90
    for item in helpitems:
        hstr = item[0].upper() + ':  ' + item[1]
        help_str += '\n'+'\n\t '.join(textwrap.wrap(hstr, basev))
    return help_str


def get_classmap():

    comp_module = modules['krangpower.components']
    classmap = {}
    for item in dir(comp_module):
        classmap[item.lower()] = getattr(comp_module, item)

    return classmap


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


def _matrix_from_json(value):

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


def load_entities(path):

    classmap = get_classmap()

    with open(path, 'r') as file:
        dik = json.load(file)

    # json entity file contain jsonized objects. This means that all lists are np.matrix.tolist representation
    # and we have to convert them back.
    for entity in dik:
        for property, value in dik[entity]['properties'].items():
            if isinstance(value, list):
                dik[entity]['properties'][property] = _matrix_from_json(value)

    dicky = {}

    for entity_name in dik:
        elcls = classmap[dik[entity_name]['type']]
        if elcls.isnamed():
            dicky[entity_name] = elcls(entity_name, xml_rep=dik[entity_name]['properties'])
        else:
            dicky[entity_name] = elcls(dik[entity_name]['properties'])

    return dicky


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
                dik[entity]['properties'][prop] = _matrix_from_json(value)

            # todo give it a unit measure

    return dik