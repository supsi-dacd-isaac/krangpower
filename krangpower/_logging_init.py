# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import logging
import os
from logging.handlers import RotatingFileHandler

from ._config_loader import GLOBAL_LOG_LEVEL, DEFAULT_ENH_NAME, MAX_LOG_SIZE_MB, MAX_HIST_LOG
from . import _config_loader as cl

_MiB = 2 ** 20


def _create_main_logger():
    logformat = '%(asctime)s - %(levelname)s (%(funcName)s) -------------> %(message)s'
    main_logger = logging.getLogger('krangpower')
    main_logger.setLevel(GLOBAL_LOG_LEVEL)
    logformatter = logging.Formatter(logformat)

    # streamhandler
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARN)
    ch.setFormatter(logformatter)
    main_logger.addHandler(ch)

    setattr(main_logger, 'official_formatter', logformatter)

    return main_logger


def _create_command_debug_logger(name):
    logformat = '%(asctime)s - %(message)s'
    cmd_logger = logging.getLogger(name)
    cmd_logger.setLevel(GLOBAL_LOG_LEVEL)
    logformatter = logging.Formatter(logformat)

    setattr(cmd_logger, 'official_formatter', logformatter)

    return cmd_logger


def _create_bare_command_logger(name):
    logformat = '%(message)s'
    bcmd_logger = logging.getLogger(name)
    bcmd_logger.setLevel(logging.DEBUG)
    logformatter = logging.Formatter(logformat)

    setattr(bcmd_logger, 'official_formatter', logformatter)

    return bcmd_logger


clog = _create_command_debug_logger(DEFAULT_ENH_NAME)
mlog = _create_main_logger()
bclog = _create_bare_command_logger('bare_' + DEFAULT_ENH_NAME)


def add_rotfilehandler(logger, path):

    logformatter = logger.official_formatter
    # this property is expected because it is set in the creation functions

    try:
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        fh = RotatingFileHandler(path, maxBytes=MAX_LOG_SIZE_MB * _MiB, backupCount=MAX_HIST_LOG)
        fh.setFormatter(logformatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)
    except PermissionError:
        # this is handled to the logger itself for warning
        logger.warning('Permission to write log file denied')


def add_comfilehandler(logger, path):

    logformatter = logger.official_formatter
    # this property is expected because it is set in the creation functions

    try:
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        fh = logging.FileHandler(path, mode='w')
        fh.setFormatter(logformatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)
    except PermissionError:
        # this is handled to the logger itself for warning
        logger.warning('Permission to write log file denied')


def remove_filehandlers(logger):
    newhandlers = [lg for lg in logger.handlers if not isinstance(lg, logging.FileHandler)]
    logger.handlers = newhandlers


def set_log_level(lvl):
    cl.GLOBAL_LOG_LEVEL = lvl
    clog.setLevel(lvl)
    mlog.setLevel(lvl)


def get_log_level():
    return cl.GLOBAL_LOG_LEVEL
