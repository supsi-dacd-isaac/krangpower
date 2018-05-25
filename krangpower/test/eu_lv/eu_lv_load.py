from krangpower import Krang, from_json

UM = Krang.get_unit_registry()

eulv = from_json("./eu_lv.json")

eulv.set(basefrequency=50.0 * UM.Hz)

eulv.set(time=[0, 0], number=1, stepsize=1.0 * UM.min)
eulv.solve()
p_1 = eulv.export('powers')
print(p_1)

eulv.set(time=[0, 60*565], number=1, stepsize=1.0 * UM.min)
eulv.solve()
p_566 = eulv.export('powers')
print(p_566)

eulv.set(time=[0, 60*1439], number=1, stepsize=1.0 * UM.min)
eulv.solve()
v_1440 = eulv.export('voltages')
print(v_1440)

eulv.set(time=[0, 0], number=1440, stepsize=1.0 * UM.min)
eulv.solve()
