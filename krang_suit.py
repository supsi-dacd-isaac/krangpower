import re
from functools import singledispatch

import networkx as nx
import pandas
from tqdm import tqdm

import components as od
from components import um
from enhancer import OpendssdirectEnhancer

_elk = 'el'


class OEKrang:

    def __init__(self):
        self.oe = OpendssdirectEnhancer()
        self.up_to_date = False
        self.last_gr = None
        self.gen_echo = False
        self.named_entities = []

        self.set(mode='duty')

    def initialize(self, name, vsource: od.Vsource):
        master_string = 'clear\nnew object = circuit.' + name + ' '
        main_source_dec = vsource.fcs(buses=('sourcebus', 'earth'))
        main_dec = re.sub('.*source ', '', main_source_dec)
        main_dec = re.sub('bus2=[^ ]+ ', '', main_dec)
        master_string += '\n~ ' + main_dec + '\n'

        self.command(master_string)
        self.command('solve')

    def __getitem__(self, item):
        # OEShell['bus.nname'] gets the component or a list of components of the bus
        # OEShell[('bus.nname',)] gets a view
        # implemented outside for single dispatching
        return _oe_getitem(item, self)

    def __getattr__(self, item):
        # OEShell.lines gets you a view of the lines
        raise NotImplementedError

    def __add__(self, other):
        assert other.isnamed()
        assert other.name != ''
        self.command(other.fcs())
        self.named_entities.append(other)
        return self

    def __bool__(self):
        return self.up_to_date

    def command(self, cmd_str: str, echo=False):
        if echo or self.gen_echo:
            print(cmd_str)
        rslt = self.oe.utils.run_command(cmd_str)
        self.up_to_date = False
        return rslt

    def set(self, **opts_vals):
        for option, value in opts_vals.items():
            self.command('set {0}={1}'.format(option, value))

    def snap(self):
        self.command('set mode=snap\nsolve\nset mode=duty')

    def _drag_solve(self):
        nmbr = self.oe.Solution.Number()
        self.oe.Solution.Number(1)
        v = pandas.DataFrame(
            columns=[x.lower() for x in self.oe.Circuit.YNodeOrder()])
        i = pandas.DataFrame(
            columns=[x.lower() for x in self.oe.Circuit.YNodeOrder()])

        for step in tqdm(range(nmbr)):
            self.command('solve')
            v = v.append(self.oe.Circuit.YNodeVArray(), ignore_index=True)
            i = i.append(self.oe.Circuit.YCurrents(), ignore_index=True)

        self.oe.Solution.Number(nmbr)

        return v, i

    @property
    def name(self):
        return self.oe.Circuit.Name()

    @property
    def graph(self):

        if self:
            return self.last_gr
        else:
            gr = nx.MultiGraph()
            ns = self.oe.Circuit.AllElementNames()
            for name in ns:
                try:
                    buses = self.oe[name].BusNames()
                except TypeError:
                    continue

                if len(buses) == 1:
                    gr.add_node(buses[0], **{_elk: self.oe[name]})
                elif len(buses) == 2:
                    gr.add_edge(buses[0], buses[1], **{_elk: self.oe[name]})
                else:
                    raise IndexError('Buses were > 2. This is a big problem.')

            self.up_to_date = True
            self.last_gr = gr

            return gr

    @property
    def bus_coords(self):
        bp = {}
        for bn in self.oe.Circuit.AllBusNames():  # todo enhance the name getting
            if self.oe['bus.' + bn].Coorddefined():
                bp[bn] = (self.oe[bn].X(), self.oe[bn].Y())
            else:
                bp[bn] = None
        return bp

    @staticmethod
    def _bus_resolve(bus_descriptor: str):

        bus_descriptor.replace('bus.', '')
        tkns = bus_descriptor.split('.')

        bus = tkns[0]
        terminals = tuple(int(x) for x in tkns[1:])

        return bus, terminals


@singledispatch
def _oe_getitem(item, oeshell):
    # no default implementation
    raise TypeError('Invalid identificator passed. You can specify fully qualified element names as str, or bus/'
                    'couples of buses as tuples of str.')


@_oe_getitem.register(str)
def _(item, oeshell):
    return oeshell.oe[item]


@_oe_getitem.register(tuple)
def _(item, oeshell):
    assert len(item) <= 2
    bustermtuples = map(oeshell._bus_resolve, item)
    return _BusView(oeshell, list(bustermtuples))


class _BusView:

    def __init__(self, oek: OEKrang, bustermtuples):
        self.btt = bustermtuples
        self.tp = dict(bustermtuples)
        self.buses = tuple(self.tp.keys())
        self.oek = oek

        buskey = 'buses'
        tkey = 'terminals'

        if len(self.buses) == 1:
            try:
                self.content = self.oek.graph.nodes[self.buses][_elk]
                # it's ok that a KeyError by both indexes is caught in the same way
            except KeyError:
                self.content = []
        elif len(self.buses) == 2:
            self.content = list(nx.get_edge_attributes(self.oek.graph.subgraph(self.buses), _elk).values())
        else:
            raise ValueError

        self.fcs_kwargs = {buskey: self.buses, tkey: self.tp}

    def __add__(self, other):
        assert not other.isnamed()
        assert other.name != ''
        self.oek.command(other.fcs(**self.fcs_kwargs))
        return _BusView(self.oek, self.btt)

    def __getitem__(self, item):
        return self.content[item]

    def __str__(self):
        return '<BusView' + str(self.buses) + ':' + str([x.Name() for x in self.content]) + '>'

    def __repr__(self):
        return self.__str__()


moe = OEKrang()
moe.gen_echo = True
pish = od.LineCode_S('thefaggy', r0=100.0 * um.ohm / um.unitlength)
posh = od.LineCode_A('thefiggy', units='m')



moe.initialize('myckt',  od.Vsource().aka('source'))
moe + pish
moe + posh
moe['sourcebus', 'a'] + od.Line(length=120).aka('theline') * posh
moe[('a', 'b')] + od.Line(units='m', length=0.12 * um.km).aka('theotherline') * posh

# print(moe['line.theotherline']['rmatrix'])
yyyy = od.Line(length=120).aka('theline') * posh

# moe['line.theotherline']['length'] = 200 * um.yard


# lish = od.CsvLoadshape(r"D:\d\ams_data\loads" + r'\Y6' + r"\node_" + str(5438) + ".csv", column_scheme={'mult': 1, 'qmult': 2}, use_actual=True, interval=15, npts=200)


pee = od.WireData('caggi', runits='m', rac=123 * um.ohm / um.m)
wee = od.WireData('aaggi')
gee = od.LineGeometry_O('faggi', nconds=2, nphases=2, x=[0, 0], h=[0, 0], units=['m', 'km']) * [pee, wee]

moe.gen_echo = True
moe + pee
moe + wee
moe + gee


print(pee.jsonize())
#
#
fofo = od.Load(kw = 1.2345 * um.MW) # * lish
moe[('b',)] + fofo.aka('wow')

tt = moe['load.wow'].unpack(verbose=False)

print(moe['line.theotherline'].unpack(verbose=True).jsonize())
print(moe['load.wow'].unpack(verbose=True).jsonize())

# moe + lish

# i, v = moe._drag_solve()
#
#
# y = moe.graph
# pass


# MOE = OpendssdirectEnhancer()
# MOE.utils.run_command("Clear")
# MOE.utils.run_command("New object = circuit.myckt bus1=sourcebus basekv=11.0 pu=1.0 angle=0.0 phases=3")
# MOE.utils.run_command("New line.line1 bus1=sourcebus bus2=miao basefreq=50.0")
# MOE.utils.run_command("New load.myload bus1=miao kw=10.0 kv=11.0 basefreq=50.0")
# u = MOE.utils.run_command("? load.myload.kw")
# MOE.utils.run_command("solve")
# print(u)
