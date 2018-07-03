# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import tqdm
from ._logging_init import set_log_level
from . import _config_loader as cl
from logging import INFO

__all__ = ['PBar']


class PBar:
    def __init__(self, iterable, level=INFO, *args, **kwargs):
        """Pbar is a modified version of the tqdm bar that accepts a logging level as an init parameter and displays the
         progress bar only if the global log level of krangpower is lower than it."""
        self.lvl = level
        self._itr = iterable
        self._args = args
        self._kwargs = kwargs
        self._tqdm = None

    def __iter__(self):
        if self.lvl >= cl.GLOBAL_LOG_LEVEL:
            self._tqdm = tqdm.tqdm(iterable=self._itr, ascii=cl.PBAR_ISASCII, *self._args, **self._kwargs)
            return self._tqdm.__iter__()
        else:
            return self._itr.__iter__()


def _main():
    original_gll = cl.GLOBAL_LOG_LEVEL

    set_log_level(10)
    print('a pbar should appear now:')
    for _ in PBar(range(10), 20, desc='appears'):
        pass

    set_log_level(20)
    print('a pbar should appear now:')
    for _ in PBar(range(10), 20, desc='appears, should be the last'):
        pass

    set_log_level(30)
    print('a pbar should NOT appear now:')
    for _ in PBar(range(10), 20, desc='SHOULD NOT APPEAR'):
        pass

    set_log_level(original_gll)


if __name__ == '__main__':
    _main()
