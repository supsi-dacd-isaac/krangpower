import argparse
import os.path

here = os.path.abspath(os.path.dirname(__file__))
y = argparse.ArgumentParser()

y.add_argument('subv')
args = y.parse_args()

with open(os.path.join(here, 'subversion.txt'), 'w') as f:
    f.write(args.subv)
