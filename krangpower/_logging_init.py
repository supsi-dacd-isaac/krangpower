import logging
import os
from logging.handlers import RotatingFileHandler

from ._config_loader import GLOBAL_LOG_LEVEL, DEFAULT_ENH_NAME, COMMAND_LOGPATH, MAIN_LOGPATH, MAX_LOG_SIZE_MB, MAX_HIST_LOG
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

    # filehandler
    try:
        if not os.path.exists(os.path.dirname(MAIN_LOGPATH)):
            os.makedirs(os.path.dirname(MAIN_LOGPATH))
        fh = RotatingFileHandler(MAIN_LOGPATH, maxBytes=MAX_LOG_SIZE_MB * _MiB, backupCount=MAX_HIST_LOG)
        fh.setFormatter(logformatter)
        fh.setLevel(logging.DEBUG)
        main_logger.addHandler(fh)
    except (PermissionError, OSError):
        # this is handled to the console stream
        main_logger.warning('Permission to write log file denied')

    return main_logger


def _create_command_logger(name):
    logformat = '%(asctime)s - %(message)s'
    cmd_logger = logging.getLogger(name)
    cmd_logger.setLevel(GLOBAL_LOG_LEVEL)
    logformatter = logging.Formatter(logformat)

    # filehandler
    try:
        if not os.path.exists(os.path.dirname(COMMAND_LOGPATH)):
            os.makedirs(os.path.dirname(COMMAND_LOGPATH))
        fh = RotatingFileHandler(COMMAND_LOGPATH, maxBytes=MAX_LOG_SIZE_MB * _MiB, backupCount=MAX_HIST_LOG)
        fh.setFormatter(logformatter)
        fh.setLevel(logging.DEBUG)
        cmd_logger.addHandler(fh)
    except PermissionError:
        # this is handled to the console stream
        cmd_logger.warning('Permission to write log file denied')

    return cmd_logger


clog = _create_command_logger(DEFAULT_ENH_NAME)
mlog = _create_main_logger()


def set_log_level(lvl):
    cl.GLOBAL_LOG_LEVEL = lvl
    clog.setLevel(lvl)
    mlog.setLevel(lvl)


def get_log_level():
    return cl.GLOBAL_LOG_LEVEL
