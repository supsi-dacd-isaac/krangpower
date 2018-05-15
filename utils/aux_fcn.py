import components as co
import json


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


def dejsonize(obj_repr: dict):

    classmap = co.get_classmap()

    # determines class
    elcls = classmap[obj_repr['type']]

    if 'path' in obj_repr.keys():
        with open(obj_repr['path'], 'r') as file:
            dik = json.load(file)
            obj_repr['properties'] = dik[obj_repr['name']]['properties']

    # restore matrices
    for prop, value in obj_repr['properties'].items():
        if isinstance(value, list):
            obj_repr['properties'][prop] = co._matricize(value)

    # add in the adjointed parameters
    obj_repr['properties'].update(obj_repr['depends'])

    # returns object
    if elcls.isnamed():
        return elcls(obj_repr['name'], xml_rep=obj_repr['properties'])
    else:
        return elcls(obj_repr['properties'])