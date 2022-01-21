def main():

    from numpy import round
    from numpy.random import random, seed
    import csv
    import os.path
    import sys

    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..\\..\\..\\krangpower'))
    import krangpower as kp

    kp.set_log_level(10)

    um = kp.UM
    src = kp.Vsource(basekv=10.0 * um.kV, pu=1.0)
    twc = kp.Krang('twc', src)
    twc['sourcebus', 'a'] << kp.Line(length=20.0 * um.m).aka('l1')

    # -------------------------------------------------------
    # 3-winding transformer
    trf = kp.Transformer(windings=3,
                         conns=['delta', 'wye', 'wye'],
                         kvs=[10.0, 2.0, 3.0] * um.kV,
                         kvas=[100.0, 20.0, 85.0] * um.kW,
                         taps=[1.0, 1.0, 1.0],
                         pctrs=[1.0, 1.0, 1.0])
    twc['a', 'b', 'c'] << trf.aka('trf')  # the first winding is considered primary
    twc['b', 'box_element'] << kp.Line(length=10.0 * um.m).aka('lbb')
    twc['c', 'cc'] << kp.Line(length=15.0 * um.m).aka('lcc')

    # -------------------------------------------------------
    # live csv loadshape creation
    with open('ls.csv', 'w') as csvfile:
        lshwriter = csv.writer(csvfile, delimiter=',')
        lshwriter.writerow(['1.0'])
        lshwriter.writerow(['1.1'])
        lshwriter.writerow(['1.2'])
        lshwriter.writerow(['0.9'])
        lshwriter.writerow(['0.2'])
        lshwriter.writerow(['1.8'])

    ls = kp.CsvLoadshape('simple_lsh', 'ls.csv', interval=2 * um.min)
    twc['cc', ] << kp.Load(kv=3.0 * um.kV, kw=35.0 * um.kW).aka('loadhi') * ls

    # -------------------------------------------------------
    # a simple fourq
    class MyDM(kp.DecisionModel):
        def decide_pq(self, oek, mynode):
            # our "smart" load
            seed(945021)  # for repeatability
            return -round(80.0 * random(), decimals=3) * um.kW,\
                   -round(10.0 * random(), decimals=3) * um.kW

    fq = kp.FourQ(kv=2.0 * um.kV)
    twc['box_element', ] << fq.aka('myfq') * MyDM()

    twc.set(number=6, stepsize=2 * um.min)

    twc.drag_solve()

    # -------------------------------------------------------
    # GraphView demonstration

    # these views report the voltages of the last solution and the base voltages
    bus_volts = kp.gv.BusVoltageView(twc)
    bus_base_volts = kp.gv.BaseVoltageView(twc, [2.0, 3.0, 10.0] * um.kV)

    # what are the buses whose voltage is < 0.98 of their nominal one?
    v_limit = 0.98
    violated_bus = {x: (bus_volts[x] / bus_base_volts[x]).to_base_units()
                    for x in bus_volts.nodes
                    if bus_volts[x] / bus_base_volts[x] < v_limit}

    print(violated_bus)

    kp.set_log_level(30)


if __name__ == '__main__':
    main()
