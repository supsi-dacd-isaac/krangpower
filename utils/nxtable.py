import networkx as nx
import csv


class NxTable:
    """A 2-index data container settable and gettable with 2-tuples of any hashable. The inner data structure is based
    on Networkx."""
    def __init__(self):
        self._data = nx.DiGraph()

    def __setitem__(self, key, value):
        self._check_key(key)
        self._data.add_edge(*key, data=value)

    def __getitem__(self, item):
        self._check_key(item)
        return self._data[item[0]][item[1]]['data']

    def from_csv(self, path):
        with open(path, 'r') as f:
            cntnt = csv.reader(f)

            for row in cntnt:
                self[row[0], row[1]] = row[2]

    @staticmethod
    def _check_key(key):
        try:
            assert isinstance(key, tuple)
            assert len(key) == 2
        except AssertionError:
            raise IndexError('Invalid key {0}'.format(key))


def _main():
    n = NxTable()
    n['foo', 2] = 3
    n['bar', 'goo'] = tuple

    print(n['bar', 'goo'])
    print(n['foo', 2])


if __name__ == '__main__':
    _main()
