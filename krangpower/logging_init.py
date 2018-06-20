import logging
import os
from logging.handlers import RotatingFileHandler

from .config_loader import _GLOBAL_LOG_LEVEL, CONFIG, _DEFAULT_NAME, _COMMAND_LOGPATH, _MAIN_LOGPATH, _MAX_LOG_SIZE
from . import config_loader as cl

MiB = 2**20


def _create_main_logger():
    logformat = '%(asctime)s - %(levelname)s (%(funcName)s) -------------> %(message)s'
    main_logger = logging.getLogger('krangpower')
    main_logger.setLevel(_GLOBAL_LOG_LEVEL)
    logformatter = logging.Formatter(logformat)

    # streamhandler
    _ch = logging.StreamHandler()
    _ch.setLevel(logging.WARN)
    _ch.setFormatter(logformatter)
    main_logger.addHandler(_ch)

    # filehandler
    try:
        if not os.path.exists(os.path.dirname(_MAIN_LOGPATH)):
            os.makedirs(os.path.dirname(_MAIN_LOGPATH))
        fh = RotatingFileHandler(_MAIN_LOGPATH, maxBytes=_MAX_LOG_SIZE*MiB, backupCount=1)
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
    cmd_logger.setLevel(_GLOBAL_LOG_LEVEL)
    logformatter = logging.Formatter(logformat)

    # filehandler
    try:
        if not os.path.exists(os.path.dirname(_COMMAND_LOGPATH)):
            os.makedirs(os.path.dirname(_COMMAND_LOGPATH))
        fh = RotatingFileHandler(_COMMAND_LOGPATH, maxBytes=_MAX_LOG_SIZE*MiB, backupCount=1)
        fh.setFormatter(logformatter)
        fh.setLevel(logging.DEBUG)
        cmd_logger.addHandler(fh)
    except PermissionError:
        # this is handled to the console stream
        cmd_logger.warning('Permission to write log file denied')

    return cmd_logger


_clog = _create_command_logger(_DEFAULT_NAME)
_mlog = _create_main_logger()


def set_log_level(lvl):
    cl._GLOBAL_LOG_LEVEL = lvl
    _clog.setLevel(lvl)
    _mlog.setLevel(lvl)


def get_log_level():
    return cl._GLOBAL_LOG_LEVEL
