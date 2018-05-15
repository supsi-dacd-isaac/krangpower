import re
from functools import singledispatch

import networkx as nx
import pandas
import json
from tqdm import tqdm

import components as co
from components import um, config
import utils.aux_fcn as au
from enhancer import OpendssdirectEnhancer
from logging.handlers import RotatingFileHandler
import logging
import os.path

_elk = 'el'


class KrangSuit:

    def __init__(self):
        self.oe = OpendssdirectEnhancer()
        self.up_to_date = False
        self.last_gr = None
        self.gen_echo = False
        self.named_entities = []
        self.ai_list = []
        self.clog = None

    def initialize(self, name, vsource: co.Vsource):
        self.clog = _create_command_logger(name)
        self.command('clear')
        master_string = 'new object = circuit.' + name + ' '
        main_source_dec = vsource.fcs(buses=('sourcebus', 'earth'))
        main_dec = re.sub('.*source ', '', main_source_dec)
        # main_dec = re.sub('bus2=[^ ]+ ', '', main_dec)
        master_string += ' ' + main_dec + '\n'

        self.command(master_string)
        self.set(mode='duty')
        self.command('makebuslist')

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

    def command(self, cmd_str: str, echo=True):
        rslt = self.oe.utils.run_command(cmd_str)
        self.up_to_date = False
        if self.gen_echo and echo:
            print(cmd_str)
            if rslt != '':
                print(str('\033[94m' + rslt + '\033[0m'))

        self.clog.debug('[' + cmd_str.replace('\n', '\n' + ' '*(30 + len(self.name)))
                        + ']-->[' + rslt.replace('\n', '') + ']')

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

        for _ in tqdm(range(nmbr)):
            for ai_el in self.ai_list:
                self.command(ai_el.element.fus(self, ai_el.name))

            self.command('solve', echo=False)
            v = v.append(self.oe.Circuit.YNodeVArray(), ignore_index=True)
            i = i.append(self.oe.Circuit.YCurrents(), ignore_index=True)

        self.oe.Solution.Number(nmbr)

        return v, i

    def make_json_dict(self):
        master_dict = {'cktname': self.name, 'elements': {}}

        for nm in self.oe.Circuit.AllElementNames():
            master_dict['elements'][nm.split('.')[1]] = self[nm].unpack().jsonize()
            master_dict['elements'][nm.split('.')[1]]['topological'] = self[nm].topological

        for ne in self.named_entities:
            master_dict['elements'][ne.name] = ne.jsonize()

        return master_dict

    def save_json(self, path):
        with open(path, 'w') as ofile:
            json.dump(self.make_json_dict(), ofile, indent=4)

    @property
    def name(self):
        return self.oe.Circuit.Name()

    @property
    def graph(self):

        if False:
            return self.last_gr
        else:
            gr = nx.MultiGraph()
            ns = self.oe.Circuit.AllElementNames()
            for name in ns:
                try:
                    buses = self.oe[name].BusNames()
                except TypeError:
                    continue

                #todo encode term perms in the graph
                if len(buses) == 1:
                    bs, _ = self._bus_resolve(buses[0])
                    gr.add_node(bs, **{_elk: self.oe[name]})
                elif len(buses) == 2:
                    bs0, _ = self._bus_resolve(buses[0])
                    bs1, _ = self._bus_resolve(buses[1])
                    gr.add_edge(bs0, bs1, **{_elk: self.oe[name]})
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


def _create_command_logger(name):
    logformat = '%(asctime)s - %(name)s - %(message)s'
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logformatter = logging.Formatter(logformat)

    # filehandler
    try:
        logpath = os.path.join(os.getenv('APPDATA'), config.get('log_file', 'commands_log_path'))
        if not os.path.exists(os.path.dirname(logpath)):
            os.makedirs(os.path.dirname(logpath))
        fh = RotatingFileHandler(logpath, maxBytes=2e6, backupCount=0)
        fh.setFormatter(logformatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)
    except PermissionError:
        # this is handled to the console stream
        logger.warning('Permission to write log file denied')

    return logger


class _DepGraph(nx.DiGraph):
    """Simple extension of nx.Digraph created to reverse-walk a dependency tree in order to declare the entities in the
    right order."""
    @property
    def leaves(self):
        return [x for x in self.nodes() if self.out_degree(x) == 0]

    def trim(self):
        self.remove_nodes_from(self.leaves)


class _BusView:

    def __init__(self, oek: KrangSuit, bustermtuples):
        self.btt = bustermtuples
        self.tp = dict(bustermtuples)
        self.buses = tuple(self.tp.keys())
        self.nb = len(self.buses)
        self.oek = oek

        buskey = 'buses'
        tkey = 'terminals'

        if len(self.buses) == 1:
            try:
                self.content = self.oek.graph.nodes[self.buses[0]][_elk]
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

        # remember ai elements
        if other.isai():
            self.oek.ai_list.append(other)

        return _BusView(self.oek, self.btt)

    def __getitem__(self, item):
        return self.content[item]

    def __getattr__(self, item):
        raise NotImplementedError

    def __str__(self):
        return '<BusView' + str(self.buses) + '>'

    def __repr__(self):
        return self.__str__()


def from_json(path):

    # load all entities
    with open(path, 'r') as ofile:
        master_dict = json.load(ofile)

    # init the krang with the source, then remove it
    l_ckt = KrangSuit()
    l_ckt.initialize(master_dict['cktname'], au.dejsonize(master_dict['elements']['source']))
    del master_dict['elements']['source']

    # reconstruction of dependency graph
    dep_graph = _DepGraph()
    for jobj in master_dict['elements'].values():
        if jobj['depends'] == {} or all([d == '' for d in jobj['depends'].values()]):
            dep_graph.add_node(jobj['name'])
        else:
            for dvalue in jobj['depends'].values():
                if isinstance(dvalue, list):
                    for dv in dvalue:
                        if dv != '':
                            dep_graph.add_edge(jobj['name'], dv)
                else:
                    if dvalue != '':
                        dep_graph.add_edge(jobj['name'], dvalue)

    while dep_graph.leaves:
        for nm in dep_graph.leaves:
            jobj = master_dict['elements'][nm]
            dssobj = au.dejsonize(jobj)
            if dssobj.isnamed():
                l_ckt + dssobj
            else:
                l_ckt[tuple(jobj['topological'])] + dssobj.aka(jobj['name'])
                # l_ckt.command(dssobj.aka(jobj['name']).fcs(buses=jobj['topological']))
        dep_graph.trim()

    return l_ckt



def _main():

    moe = KrangSuit()
    moe.gen_echo = True
    pish = co.LineCode_S('thefaggy', r0=100.0 * um.ohm / um.unitlength)
    posh = co.LineCode_A('thefiggy', units='m')

    moe.initialize('myckt', co.Vsource(basekv=11.0 * um.kV).aka('source'))
    moe + pish
    moe + posh
    moe['sourcebus', 'a'] + co.Line(length=120 * um.unitlength).aka('theline') * posh
    moe[('a.1.3.2', 'b.3.2.1')] + co.Line(units='m', length=0.72 * um.km).aka('theotherline')

    # print(moe['line.theotherline']['rmatrix'])
    yyyy = co.Line(length=120).aka('theline') * posh

    # moe['line.theotherline']['length'] = 200 * um.yard


    lish = co.CsvLoadshape(r"D:\d\ams_data\loads" + r'\Y6' + r"\node_" + str(5438) + ".csv", column_scheme={'mult': 1, 'qmult': 2}, use_actual=True, interval=15 * um.min)


    pee = co.WireData('caggi', runits='m', rac=123 * um.ohm / um.m)
    wee = co.WireData('aaggi')
    gee = co.LineGeometry_O('faggi', nconds=2, nphases=2, x=[0, 0], h=[0, 0], units=['m', 'km']) * [pee, wee]


    moe + pee
    moe + wee
    moe + gee

    fofo = co.Load(kw=1.2345 * um.MW)  # * lish
    moe[('b',)] + fofo.aka('wow') * lish

    tt = moe['load.wow'].unpack(verbose=False)

    moe.command('makebuslist')
    print(moe.oe.Circuit.AllBusNames())
    print(moe[('b',)])
    print(moe.graph.edges)

    i,v = moe._drag_solve()

    moe.save_json(r'D:\temp\krang.json')
    # cs = from_json(r'D:\temp\krang.json')



if __name__ == '__main__':
    _main()
