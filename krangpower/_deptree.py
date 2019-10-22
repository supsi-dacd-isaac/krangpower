# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import networkx as nx

__all__ = ['DepTree']


class DepTree(nx.DiGraph):
    """Simple extension of nx.Digraph created to reverse-walk a dependency branching in order to declare the entities
    in the right order."""

    @property
    def leaves(self):
        return [x for x in self.nodes() if not nx.descendants(self, x)]

    def trim(self):
        self.remove_nodes_from(self.leaves)

    def recursive_prune(self):

        # dependency non-circularity check
        cy = list(nx.simple_cycles(self))
        try:
            assert len(cy) == 0
        except AssertionError:
            print('Found circular dependency(s) in the graph!')
            for cycle in cy:
                print(cycle)
            raise ValueError('Cannot recursively prune a directed graph with cycles')

        # actual prune generator
        while self.leaves:
            yield self.leaves
            self.trim()
