import re
from functools import singledispatch

import networkx as nx
import pandas
import json
from tqdm import tqdm
import copy
import numpy as np

import components as co
from components import um, config, _pint_qty_type
import utils.aux_fcn as au
from enhancer import OpendssdirectEnhancer
from logging.handlers import RotatingFileHandler
import logging
import os.path
import busquery as bq

_elk = 'el'


class Krang:

    def __init__(self):
        self.oe = OpendssdirectEnhancer()
        self.up_to_date = False
        self.last_gr = None
        self.named_entities = []
        self.ai_list = []
        self.clog = None
        self.com = ''

    def initialize(self, name, vsource: co.Vsource):
        self.clog = _create_command_logger(name)
        self.command('clear')
        master_string = 'new object = circuit.' + name + ' '
        vsource.name = 'source'  # we force the name, so an already named vsource can be used as source.
        main_source_dec = vsource.fcs(buses=('sourcebus', 'sourcebus.0.0.0'))
        main_dec = re.sub('New vsource\.source ', '', main_source_dec)
        main_dec = re.sub('bus2=[^ ]+ ', '', main_dec)
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
        """Krang.item, aside from retrieving the built.in attributes, wraps by default the calls to opendssdirect's
        'class_to_dataframe' utility function. These are accessible via capital letter calls. Both singular and plural
        are accepted (e.g., 'Line' or 'Lines', 'RegControl' , 'Regcontrol', 'RegControls', but not 'transformers')"""
        try:
            assert item[0].isupper()
            dep_item = re.sub('s$', '', item).lower()
            return self.oe.utils.class_to_dataframe(dep_item)
        except (AssertionError, NotImplementedError):
            raise AttributeError('{0} is neither a valid attribute nor a valid identifier for the class-views.'.format(
                item
            ))

    def __lshift__(self, other):
        try:
            assert other.isnamed()
            self.named_entities.append(other)
        except AssertionError:
            assert other.isabove()

        assert other.name != ''
        self.command(other.fcs())
        return self

    def __bool__(self):
        return self.up_to_date

    def command(self, cmd_str: str, echo=True):
        rslt = self.oe.utils.run_command(cmd_str)
        self.up_to_date = False
        if echo:
            self.com += cmd_str + '\n'
            self.clog.debug('[' + cmd_str.replace('\n', '\n' + ' '*(30 + len(self.name)))
                            + ']-->[' + rslt.replace('\n', '') + ']')

        return rslt

    def set(self, **opts_vals):
        for option, value in opts_vals.items():
            if isinstance(value, _pint_qty_type):
                vl = value.to(um.parse_units(co.default_settings['units'][option])).magnitude
            else:
                vl = value
            self.command('set {0}={1}'.format(option, vl))

    def get(self, *opts):

        assert all([x in list(co.default_settings['values'].keys()) for x in opts])
        r_opts = {opt: self.command('get {0}'.format(opt), echo=False).split('!')[0] for opt in opts}

        for op in r_opts.keys():
            tt = type(co.default_settings['values'][op])
            if tt is list:
                r_opts[op] = eval(r_opts[op])  # lists are literals like [13]
            else:
                r_opts[op] = tt(r_opts[op])
            if op in co.default_settings['units'].keys():
                r_opts[op] *= um.parse_units(co.default_settings['units'][op])

        return r_opts

    def snap(self):
        self.command('set mode=snap\nsolve\nset mode=duty')

    def drag_solve(self):
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

    def solve(self):
        self.command('solve')

    def make_json_dict(self):
        master_dict = {'cktname': self.name, 'elements': {}, 'settings': {}}

        for Nm in self.oe.Circuit.AllElementNames():
            nm = Nm.lower()
            master_dict['elements'][nm] = self[nm].unpack().jsonize()
            master_dict['elements'][nm]['topological'] = self[nm].topological

        for ne in self.named_entities:
            master_dict['elements'][ne.fullname] = ne.jsonize()

        # options
        opts = self.get(*list(co.default_settings['values'].keys()))
        for on, ov in opts.items():
            if isinstance(ov, _pint_qty_type):
                opts[on] = ov.to(um.parse_units(co.default_settings['units'][on])).magnitude
                if isinstance(opts[on], (np.ndarray, np.matrix)):
                    opts[on] = opts[on].tolist()

        master_dict['settings']['values'] = opts
        master_dict['settings']['units'] = co.default_settings['units']

        return master_dict

    def save_json(self, path):
        with open(path, 'w') as ofile:
            json.dump(self.make_json_dict(), ofile, indent=4)

    def save_dss(self, path):
        with open(path, 'w') as ofile:
            ofile.write(self.com)

    @property
    def name(self):
        return self.oe.Circuit.Name()

    @property
    def graph(self):

        def _update_node(self, gr, bs, name):
            try:
                exel = gr.nodes[bs][_elk]
            except KeyError:
                gr.add_node(bs, **{_elk: [self.oe[name]]})
                return
            exel.append(self.oe[name])
            return

        def _update_edge(self, gr, ed, name):
            try:
                exel = gr.edges[ed][_elk]
            except KeyError:
                gr.add_edge(*ed, **{_elk: [self.oe[name]]})
                return
            exel.append(self.oe[name])
            return

        if False:
            return self.last_gr
        else:
            gr = nx.Graph()
            ns = self.oe.Circuit.AllElementNames()
            for name in ns:
                try:
                    buses = self.oe[name].BusNames()
                except TypeError:
                    continue

                #todo encode term perms in the graph
                if len(buses) == 1:
                    bs, _ = self._bus_resolve(buses[0])

                    _update_node(self, gr, bs, name)

                    # gr.add_node(bs, **{_elk: self.oe[name]})
                elif len(buses) == 2:
                    bs0, _ = self._bus_resolve(buses[0])
                    bs1, _ = self._bus_resolve(buses[1])

                    _update_edge(self, gr, (bs0, bs1), name)

                    # gr.add_edge(bs0, bs1, **{_elk: self.oe[name]})
                else:
                    raise IndexError('Buses were > 2. This is a big problem.')

            self.up_to_date = True
            self.last_gr = gr
        return gr

    @property
    def bus_coords(self):
        bp = {}
        for bn in self.oe.Circuit.AllBusNames():
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

    def __init__(self, oek: Krang, bustermtuples):
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

    def __lshift__(self, other):
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
        try:
            # attributes requested via getattr are searched in busquery
            f = bq.get_fun(item)
        except KeyError:
            raise AttributeError('Attribute/query function {0} is not implemented')

        if self.nb == 1:
            return f(self.oek, self, self.buses[0])
        elif self.nb == 2:
            return f(self.oek, self, self.buses)
        else:
            raise AttributeError

    def __str__(self):
        return '<BusView' + str(self.buses) + '>'

    def __repr__(self):
        return self.__str__()


def from_json(path):

    # load all entities
    with open(path, 'r') as ofile:
        master_dict = json.load(ofile)

    # init the krang with the source, then remove it
    l_ckt = Krang()
    l_ckt.initialize(master_dict['cktname'], au.dejsonize(master_dict['elements']['vsource.source']))
    del master_dict['elements']['vsource.source']

    # load and declare options
    opt_dict = master_dict['settings']
    for on, ov in opt_dict['values'].items():
        # we try in any way to see if the value is the same as the default and, if so, we continue
        # todo there is an edge case where the value of a measured quantity is the same, but the unit is different
        if ov == co.default_settings['values'][on]:
            continue
        try:
            if np.isclose(ov, co.default_settings['values'][on]):
                continue
        except ValueError:
            if np.isclose(ov, co.default_settings['values'][on]).all():
                continue
        except TypeError:
            pass

        try:
            if ov.lower() == co.default_settings['values'][on].lower():
                continue
        except AttributeError:
            pass

        if on in opt_dict['units'].keys():
            d_ov = ov * um.parse_units(opt_dict['units'][on])
        else:
            d_ov = ov

        l_ckt.set(**{on: d_ov})

    # reconstruction of dependency graph and declarations
    dep_graph = _DepGraph()
    for jobj in master_dict['elements'].values():
        # if the element has no dependencies, we just add a node with iths name
        if jobj['depends'] == {} or all([d == '' for d in jobj['depends'].values()]):
            dep_graph.add_node(jobj['type'] + '.' + jobj['name'])
        else:
            # if an element parameter depends on another name, or a list of other names, we create all the edges
            # necessary
            for dvalue in jobj['depends'].values():
                if isinstance(dvalue, list):
                    for dv in dvalue:
                        if dv != '':
                            dep_graph.add_edge(jobj['type'] + '.' +jobj['name'], dv)
                else:
                    if dvalue != '':
                        dep_graph.add_edge(jobj['type'] + '.' + jobj['name'], dvalue)

    # we cyclically consider all "leaves", add the objects at the leaves, then trim the leaves and go on with
    # the new leaves.
    # In this way we are sure that, whenever a name is mentioned in a fcs, its entity was already declared.
    while dep_graph.leaves:
        for nm in dep_graph.leaves:
            try:
                jobj = copy.deepcopy(master_dict['elements'][nm])
            except KeyError:
                mdmod = {k.split('.')[1]: v for k, v in master_dict['elements'].items()}
                jobj = copy.deepcopy(mdmod[nm])
            dssobj = au.dejsonize(jobj)
            if dssobj.isnamed():
                l_ckt << dssobj
            elif dssobj.isabove():
                l_ckt << dssobj.aka(jobj['name'])
            else:
                l_ckt[tuple(jobj['topological'])] << dssobj.aka(jobj['name'])
                # l_ckt.command(dssobj.aka(jobj['name']).fcs(buses=jobj['topological']))
        dep_graph.trim()

    return l_ckt


def _main():

    moe = Krang()

    moe.initialize('myckt', co.Vsource(basekv=15.0 * um.kV).aka('source'))
    moe.set(number=180, stepsize='15')

    moe['sourcebus', 'a'] << co.Line(length=120 * um.unitlength).aka('theline')

    lili = co.Line(units='m', length=0.72 * um.km).aka('theotherline')

    moe << co.Monitor().aka('themon') * lili

    moe[('a.1.3.2', 'b.3.2.1')] << lili

    vls = np.matrix([15.0, 7.0]) * um.kV

    ui = co.Transformer(windings=2, kvs=vls).aka('trans')

    # moe[('b', 'c')] + ui
    #  EEEEEEEEEEEEEEEEEEEEEEEEEE

    # print(moe['line.theotherline']['rmatrix'])

    # moe['line.theotherline']['length'] = 200 * um.yard

    lish = co.CsvLoadshape('mylsh', r"D:\d\ams_data\loads" + r'\Y6' + r"\node_" + str(5413) + ".csv", column_scheme={'mult': 1, 'qmult': 2}, use_actual=True, interval=15 * um.min)

    pee = co.WireData('caggi', runits='m', rac=3 * um.ohm / um.m)
    wee = co.WireData('aaggi')
    gee = co.LineGeometry_O('faggi', nconds=2, nphases=2, x=[0, 0], h=[0.1, 0], units=['m', 'km']) * [pee, wee]



    print(gee['units'])

    moe << lish
    moe << pee
    moe << wee
    moe << gee

    au.dejsonize(lish.jsonize())

    fofo = co.Load(kw=18.2345 * um.kW, kv=15 * um.kV)  # * lish
    moe[('b',)] << fofo.aka('wow') * lish
    moe[('b',)] << fofo.aka('bow') * lish


    tt = moe['load.wow'].unpack(verbose=False)

    moe.command('makebuslist')

    i, v = moe.drag_solve()
    print(i)
    print(moe[('b',)].voltage)

    moe.save_json(r'D:\temp\krang.json')
    moe.save_dss(r'D:\temp\krang.dss')
    print(moe['vsource.source']['isc3'])
    cs = from_json(r'D:\temp\krang.json')
    # print(moe.graph.nodes)
    # print(cs.graph.nodes)


    i, v = cs.drag_solve()


if __name__ == '__main__':
    _main()
