import tqdm
from .logging_init import set_log_level
from . import config_loader as cl
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
        if self.lvl >= cl._GLOBAL_LOG_LEVEL:
            self._tqdm = tqdm.tqdm(iterable=self._itr, ascii=cl._PBAR_ISASCII, *self._args, **self._kwargs)
            return self._tqdm.__iter__()
        else:
            return self._itr.__iter__()


def main():
    original_gll = cl._GLOBAL_LOG_LEVEL

    set_log_level(10)
    for _ in PBar(range(10), 20, desc='appears'):
        pass

    set_log_level(20)
    for _ in PBar(range(10), 20, desc='appears'):
        pass

    set_log_level(30)
    for _ in PBar(range(10), 20, desc='notappears'):
        pass

    set_log_level(original_gll)


if __name__ == '__main__':
    main()
