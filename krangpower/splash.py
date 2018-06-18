import os.path
from .config_loader import krang_directory


def splash():

    with open(os.path.join(krang_directory, 'defaults/splash.txt'), 'r') as splash_file:
        print(splash_file.read())
