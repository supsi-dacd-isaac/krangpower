# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

from ._nxtable import _main as nxm
from ._pbar import _main as pbm

# this module is for inserting small edge case tests in order to increase coverage.


def do_edge_tests():
    nxm()
    pbm()
