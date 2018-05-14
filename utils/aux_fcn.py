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