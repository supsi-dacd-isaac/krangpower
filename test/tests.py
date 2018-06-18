import runpy
import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..\\'))
print(sys.path)
import krangpower as kp
kp.splash()

print('Commencing tests...')

k = os.path.dirname(os.path.realpath(__file__))
os.chdir(k)
print('Usage1...')
runpy.run_path('usage_1/us1.py', run_name='__main__')
os.chdir(k)
print('EU_lv full test...')
runpy.run_path('eu_lv/eu_lv_fulltest.py', run_name='__main__')
