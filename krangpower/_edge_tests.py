from .nxtable import _main as nxm
from .pbar import _main as pbm

# this module is for inserting small edge case tests in order to increase coverage.


def do_edge_tests():
    nxm()
    pbm()
