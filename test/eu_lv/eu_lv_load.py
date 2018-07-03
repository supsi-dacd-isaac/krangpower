import krangpower as kp

kp.set_log_level(0)
loaded_eulv = kp.open_ckt('./eu_lv.zip')

loaded_eulv.snap()
print(loaded_eulv['bus.2'].Voltages())

