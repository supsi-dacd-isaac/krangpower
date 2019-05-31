# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import hashlib
import io
import pickle
from abc import abstractmethod

from ._config_loader import UM
from ._components import Generator, _odssrep, FcsAble, FusAble
from ._aux_fcn import termrep

__all__ = ['DecisionModel', 'FourQ']


class DecisionModel:
    eltype = 'decisionmodel'

    @abstractmethod
    def decide_pq(self, oek, mynode):
        """Takes a graph and the node where the decision model has to interpret. Returns P(active power) and Q
        (reactive power), in that order, as a tuple. P and Q are not constrained in any way one to the other."""
        pass


class FourQ(Generator, FusAble):

    def __init__(self, **parameters):
        self._dm = None
        super().__init__(**parameters)

    def update_pq(self, oek, mybus):
        assert self._dm is not None
        p, q = self._dm.decide_pq(oek, mybus)
        return p, q

    def fcs(self, **hookup):

        buses = hookup['buses']
        cname = self.name
        term_perm = hookup.get('terminals', None)

        s1 = 'New ' + 'generator' + '.' + cname  # _ splitting to allow name personalization outside dss
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
        return s1 + s3 + s2

    def fus(self, oek, myname):
        mybus = oek[myname].topological[0]
        p, q = self.update_pq(oek, mybus)
        s = 'edit generator.' + myname + ' kw=' + str(p.to(UM.kW).magnitude) + ' kvar=' + str(q.to(UM.kVA).magnitude)
        return s

    def _calchash(self):
        hash_bio = io.BytesIO()
        pickle.dump(self._dm, hash_bio, protocol=pickle.HIGHEST_PROTOCOL)
        return hashlib.md5(hash_bio.getvalue()).hexdigest()

    def jsonize(self, all_params=False, flatten_mtx=True):
        md = super().jsonize(all_params, flatten_mtx)
        md['hash'] = self._calchash()

        return md

    def define_dm(self, dm):
        assert isinstance(dm, DecisionModel)
        self._dm = dm
