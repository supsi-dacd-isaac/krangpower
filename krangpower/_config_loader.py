import configparser as cfp
import json
import logging
import os.path
import re
import platform

import pint

from ._aux_fcn import load_dictionary_json

_THISDIR = os.path.dirname(os.path.realpath(__file__))
_WIN_ENV_VAR_REGEX = re.compile('(%)([^%]+)(%)')
_LINUX_ENV_VAR_REGEX = re.compile('(\$)([^(/|\\)]+)')

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
INTERFACE_METHODS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'interfaces'))
UNIT_MEASUREMENT_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'measurement_units'))
MANDATORY_UNITS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'mandatory_units'))
TREATMENTS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'treatments'))
INTERF_SELECTORS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'interface_selectors'))
DEFAULT_SETTINGS_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'default_settings'))
DEFAULT_ENTITIES_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'default_entities'))
ASSOCIATION_TYPES_PATH = os.path.join(_THISDIR, CONFIG.get('data_files', 'association_types'))

GLOBAL_LOG_LEVEL = getattr(logging, CONFIG.get('log_settings', 'default_logging_level'))
MAX_LOG_SIZE_MB = CONFIG.getfloat('log_settings', 'max_log_size_mb')
ELK = CONFIG.get('misc_settings', 'graph_element_tag')
DEFAULT_KRANG_NAME = CONFIG.get('misc_settings', 'default_krang_name')
CMD_LOG_NEWLINE_LEN = CONFIG.getint('misc_settings', 'newline_cmdlog_length')
DEFAULT_ENH_NAME = CONFIG.get('misc_settings', 'default_enhancer_name')
GLOBAL_PRECISION = CONFIG.getint('precision', 'global_precision')
LSH_ZIP_NAME = CONFIG.get('misc_settings', 'inner_loadshape_zip_filename')
PBAR_ISASCII = CONFIG.getboolean('misc_settings', 'ascii_pbar')


# -------------------------------------------------------------
#  LOG PATH
# -------------------------------------------------------------
def replace_env(match):
    return os.getenv(match.group(2))


# general log
if platform.system() == 'Windows':

    basepath = CONFIG.get('log_settings', 'win_log_folder')
    basepath = _WIN_ENV_VAR_REGEX.sub(replace_env, basepath)

    MAIN_LOGPATH = os.path.join(basepath,
                                CONFIG.get('log_settings', 'general_log_name'))
elif platform.system() == 'Linux':

    basepath = CONFIG.get('log_settings', 'linux_log_folder')
    basepath = _LINUX_ENV_VAR_REGEX.sub(replace_env, basepath)

    MAIN_LOGPATH = os.path.join(basepath,
                                CONFIG.get('log_settings', 'general_log_name'))
else:
    raise OSError('Could not find a valid log path.')

# command_log
if platform.system() == 'Windows':

    basepath = CONFIG.get('log_settings', 'win_log_folder')
    basepath = _WIN_ENV_VAR_REGEX.sub(replace_env, basepath)
    # basepath = re.sub('(%)([^%]+)(%)', replace_env, basepath)

    COMMAND_LOGPATH = os.path.join(basepath,
                                   CONFIG.get('log_settings', 'commands_log_name'))
elif platform.system() == 'Linux':

    basepath = CONFIG.get('log_settings', 'linux_log_folder')
    basepath = _LINUX_ENV_VAR_REGEX.sub(replace_env, basepath)
    # basepath = re.sub('(\$)([^/|\\]+)', replace_env, basepath)

    COMMAND_LOGPATH = os.path.join(basepath,
                                   CONFIG.get('log_settings', 'commands_log_name'))
else:
    raise OSError('Could not find a valid log path.')


# -------------------------------------------------------------
#  TEMPORARY FILES PATH
# -------------------------------------------------------------
if platform.system() == 'Windows':
    TMP_PATH = os.path.join(os.getenv('TEMP'), CONFIG.get('temp_files', 'temp_subfolder'))
elif platform.system() == 'Linux':
    TMP_PATH = os.path.join('/var/tmp', CONFIG.get('temp_files', 'temp_subfolder'))
else:
    raise OSError('Could not find a valid temp path.')

COORDS_FILE_PATH = os.path.join(TMP_PATH, 'bus_coords_fromkml.csv')

# -------------------------------------------------------------
#  UNIT MEASURE REGISTRY
# -------------------------------------------------------------
UM = pint.UnitRegistry()
UM.define('percent = 0.01 * dimensionless = pct')
UM.define('none = [generic_length] = unitlength')  # when lengths are set as none, this creates a common basis
UM.define('mt = meter')
PINT_QTY_TYPE = type(1 * UM.m)


# -------------------------------------------------------------
#  DICTIONARY OF MANDATORY UNITS
# -------------------------------------------------------------
with open(MANDATORY_UNITS_PATH, 'r') as mf:
    MANDATORY_UNITS = json.load(mf)

# -------------------------------------------------------------
#  DEFAULTS DICTIONARIES
# -------------------------------------------------------------
DEFAULT_COMP = load_dictionary_json(DEFAULT_ENTITIES_PATH)
with open(DEFAULT_SETTINGS_PATH, 'r') as f:
    DEFAULT_SETTINGS = json.load(f)

krang_directory = _THISDIR  # we expose our path, because you could want to find the test cases, etc
