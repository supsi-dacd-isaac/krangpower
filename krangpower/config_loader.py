import configparser as cfp
import json
import logging
import os.path
import platform

import pint

from krangpower.aux_fcn import load_dictionary_json

__all__ = ['UM']
_THISDIR = os.path.dirname(os.path.realpath(__file__))

# -------------------------------------------------------------
# CONFIG LOAD
# -------------------------------------------------------------
CONFIG = cfp.ConfigParser()
CONFIG.read(os.path.join(_THISDIR, 'config/krang_config.cfg'))

# -------------------------------------------------------------
# HELP LOAD
# -------------------------------------------------------------
DSSHELP = cfp.RawConfigParser()
DSSHELP.read(os.path.join(_THISDIR, 'config/DSSHelp.cfg'))


# -------------------------------------------------------------
# CONSTANTS
# -------------------------------------------------------------
_INTERFACE_METHODS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'interfaces'))
_UNIT_MEASUREMENT_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'measurement_units'))
_TREATMENTS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'treatments'))
_INTERF_SELECTORS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'interface_selectors'))
_DEFAULT_SETTINGS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'default_settings'))
_DEFAULT_ENTITIES_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'default_entities'))
_ASSOCIATION_TYPES_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'association_types'))

_GLOBAL_LOG_LEVEL = getattr(logging, CONFIG.get('misc_settings', 'default_logging_level'))
_ELK = CONFIG.get('misc_settings', 'graph_element_tag')
_DEFAULT_KRANG_NAME = CONFIG.get('misc_settings', 'default_krang_name')
_CMD_LOG_NEWLINE_LEN = CONFIG.getint('misc_settings', 'newline_cmdlog_length')
_DEFAULT_NAME = CONFIG.get('misc_settings', 'default_enhancer_name')


# -------------------------------------------------------------
#  LOG PATH
# -------------------------------------------------------------
# general log
if platform.system() == 'Windows':
    _MAIN_LOGPATH = os.path.join(os.getenv('APPDATA'), CONFIG.get('log_file', 'log_folder'),
                                 CONFIG.get('log_file', 'general_log_path'))
elif platform.system() == 'Linux':
    _MAIN_LOGPATH = os.path.join('/var/log', CONFIG.get('log_file', 'log_folder'),
                                 CONFIG.get('log_file', 'general_log_path'))
else:
    raise OSError('Could not find a valid log path.')

# command_log
if platform.system() == 'Windows':
    _COMMAND_LOGPATH = os.path.join(os.getenv('APPDATA'), CONFIG.get('log_file', 'log_folder'),
                                    CONFIG.get('log_file', 'commands_log_path'))
elif platform.system() == 'Linux':
    _COMMAND_LOGPATH = os.path.join('/var/log', CONFIG.get('log_file', 'log_folder'),
                                    CONFIG.get('log_file', 'commands_log_path'))
else:
    raise OSError('Could not find a valid log path.')


# -------------------------------------------------------------
#  TEMPORARY FILES PATH
# -------------------------------------------------------------
if platform.system() == 'Windows':
    _TMP_PATH = os.path.join(os.getenv('TEMP'), CONFIG.get('temp_folder', 'temp_folder'))
elif platform.system() == 'Linux':
    _TMP_PATH = os.path.join('/var/tmp', CONFIG.get('temp_folder', 'temp_folder'))
else:
    raise OSError('Could not find a valid temp path.')


# -------------------------------------------------------------
#  UNIT MEASURE REGISTRY
# -------------------------------------------------------------
UM = pint.UnitRegistry()
UM.define('percent = 0.01 * dimensionless = pct')
UM.define('none = [generic_length] = unitlength')  # when lengths are set as none, this creates a common basis
UM.define('mt = meter')
_PINT_QTY_TYPE = type(1 * UM.m)


# -------------------------------------------------------------
#  DEFAULTS DICTIONARIES
# -------------------------------------------------------------
DEFAULT_COMP = load_dictionary_json(_DEFAULT_ENTITIES_PATH)
with open(_DEFAULT_SETTINGS_PATH, 'r') as f:
    DEFAULT_SETTINGS = json.load(f)
