# ,---------------------------------------------------------------------------,
# |  This module is part of the krangpower electrical distribution simulation |
# |  suit by Federico Rosato <federico.rosato@supsi.ch> et al.                |
# |  Please refer to the license file published together with this code.      |
# |  All rights not explicitly granted by the license are reserved.           |
# '---------------------------------------------------------------------------'

import os.path
from ._config_loader import krang_directory


def splash():

    with open(os.path.join(krang_directory, 'defaults/splash.txt'), 'r') as splash_file:
        print(splash_file.read())
