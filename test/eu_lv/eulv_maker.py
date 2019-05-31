import csv
import io
import os
import sys
import zipfile

import requests

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..\\..\\..\\krangpower'))
import krangpower as kp

um = kp.UM
kp.set_log_level(10)
test_dir = os.path.join(kp.TMP_PATH, 'eulvtest')


def download_extract_zip(url):
    """
    Download a ZIP file and extract its contents in memory
    yields (filename, file-like object) pairs
    """

    response = requests.get(url)
    path = os.path.join(test_dir, 'eulv_originals')
    with zipfile.ZipFile(io.BytesIO(response.content)) as thezip:
        thezip.extractall(path)

    return path

# this script loads the files for the European Low Voltage Feeder, exactly as downloadable as of 21.05.2018 from here
# http://sites.ieee.org/pes-testfeeders/resources/
# and then executes the solutions.


def make_circuit(tname, save_json, save_ckt, save_effnodes, insert_loads, simplify, refine, loadmodel='2'):

    print('Downloading and extracting data from sites.ieee.org.....')
    data_url = 'http://sites.ieee.org/pes-testfeeders/files/2017/08/European_LV_Test_Feeder_v2.zip'
    eulv_root = os.path.join(download_extract_zip(data_url), 'European_LV_CSV')
    print('Finished')

    print('Extracting load profiles...')
    lp_folder = os.path.join(eulv_root, 'Load Profiles')
    # loading loadshapes in dict
    files = os.listdir(lp_folder)
    lp_dict = dict()
    for lp in files:
        nlp = int(lp.split('_')[-1].split('.')[0])
        lp_dict[nlp] = kp.CsvLoadshape('shape_' + str(nlp), os.path.join(lp_folder, lp), {'mult': 2}, 1 * um.min, False)

    # loading linecodes in dict
    print('Extracting linecodes...')
    lc_dict = dict()
    with open(os.path.join(eulv_root, 'LineCodes.csv'), 'r') as lcf:
        rd = csv.reader(lcf)
        next(rd)  # first line is a comment
        hd = next(rd)[1:]  # header without name
        while True:
            try:
                lcs = next(rd)
                for idx, item in enumerate(lcs):
                    try:
                        lcs[idx] = int(item)
                    except ValueError:
                        try:
                            lcs[idx] = float(item)
                        except ValueError:
                            continue

                unimis = um.ohm / um.parse_units(lcs[-1])
                unicap = um.nF / um.parse_units(lcs[-1])
                for char in range(2, 6):
                    lcs[char] = lcs[char] * unimis
                for char in range(6, 8):
                    lcs[char] = lcs[char] * unicap

                lc_dict[lcs[0].replace('.', '')] = dict(zip(hd, lcs[1:-1]))  # kp.linecode(dict....)
            except StopIteration:
                break

    # loading source's kw as dict
    print('Extracting source...')
    skw = dict()
    with open(os.path.join(eulv_root, 'Source.csv'), 'r') as sf:
        for nl, content in enumerate(sf.readlines()):
            if nl in (0, 1):
                continue
            else:
                ct = content.split('=')
                skw[ct[0]] = um.parse_expression(ct[1].strip('/n'))
    skw['basekv'] = skw.pop('Voltage')
    skw['frequency'] = 50.0 * um.Hz

    # the transformer file is a nuisance
    print('Instantiating xformer...')
    TR1 = kp.Transformer(phases=3,
                         kvs=[11.0, 0.416] * um.kV,
                         kvas=[0.8, 0.8] * um.MVA,
                         conns=['delta', 'wye'],
                         xhl=4.0,
                         pctrs=[0.4] * 2,
                         sub='y')

    # lines
    print('Extracting lines profiles...')
    lines_dict = dict()
    with open(os.path.join(eulv_root, 'Lines.csv'), 'r') as lnf:
        rd = csv.reader(lnf)
        next(rd)  # first line is a comment
        hd = next(rd)  # header without name
        while True:
            try:
                rl = next(rd)
                ll = dict(zip(hd[1:], rl[1:]))
                lines_dict[rl[0]] = {}
                lines_dict[rl[0]]['kwargs'] = {}
                lines_dict[rl[0]]['Buses'] = [ll['Bus1'], ll['Bus2']]
                lines_dict[rl[0]]['kwargs']['units'] = ll['Units']
                lines_dict[rl[0]]['kwargs']['length'] = float(ll['Length']) * um.m
                lines_dict[rl[0]]['kwargs']['phases'] = len(ll['Phases'])
                lines_dict[rl[0]]['LineCode'] = ll['LineCode'].replace('.', '')
            except StopIteration:
                break

    # the data is contained in:
    # skw, loads_dict, lines_dict, lc_dict, lp_dict, transformer
    print('Creating circuit...')
    SRC = kp.Vsource(**skw)
    eulv = kp.Krang(tname, SRC)
    eulv.set(basefreq=50.0 * um.Hz, voltagebases=[11.0, 0.416] * um.kV)
    eulv.command('calcvoltagebases')

    print('Adding transformer...')
    eulv['sourcebus', '1'] << TR1.aka('TR1')

    print('Adding linecodes...')
    lc_el = {}
    for lcname, lcdata in lc_dict.items():
        lc_el[lcname] = kp.LineCode(lcname, **lcdata)
        eulv << lc_el[lcname]

    # <editor-fold desc="Customization">
    def degree_two_nodes(grph):
        return [x for x in tg.nodes if grph.degree(x) == 2]

    def degree_one_nodes(grph):
        return [x for x in tg.nodes if grph.degree(x) == 1]

    def edges_analogous(graph, e1, e2):
        if graph.edges[e1]['linecode'] != graph.edges[e2]['linecode']:
            return False
        if graph.edges[e1]['phases'] != graph.edges[e2]['phases']:
            return False
        if graph.edges[e1]['units'] != graph.edges[e2]['units']:
            return False

        return True

    if simplify or refine:
        print('Lines postelaboration...')
        import networkx as nx
        tg = nx.Graph()
        for lname, line in lines_dict.items():
            tg.add_edge(*line['Buses'],
                        linecode=line['LineCode'],
                        phases=line['kwargs']['phases'],
                        length=line['kwargs']['length'],
                        units=line['kwargs']['units'])

    if simplify:
        while True:
            exit1 = False
            tn = iter(degree_two_nodes(tg))
            while True:
                try:
                    node = next(tn)
                except StopIteration:
                    exit1 = True
                    break
                e1, e2 = tg.edges(node)
                n1, n2 = tg.neighbors(node)
                if edges_analogous(tg, e1, e2):
                    tg.add_edge(n1,
                                n2,
                                phases=tg.edges[e1]['phases'],
                                linecode=tg.edges[e1]['linecode'],
                                units=tg.edges[e1]['units'],
                                length=tg.edges[e1]['length'] + tg.edges[e2]['length'])
                    tg.remove_node(node)
                    break
                else:
                    continue

            if exit1:
                break

    # further elaboration: fake-test nodes are removed
    # reading csv file
    if refine:
        with open(os.path.join(eulv_root, 'Buscoords.csv'), 'r') as csvfile:
            # creating a csv reader object
            csvreader = csv.reader(csvfile)
            # extracting each data row one by one
            next(csvreader)
            next(csvreader)

            for row in csvreader:
                bus, x, y = row
                try:
                    tg.nodes[bus]['pos'] = (float(x), float(y))
                except KeyError:
                    continue

        terminals = degree_one_nodes(tg)
        terminals_pos = {t: tg.nodes[t]['pos'] for t in terminals}
        count = 0
        for n in terminals:
            positions = nx.get_node_attributes(tg, 'pos')
            duplicate_nodes = {x: v for x, v in positions.items() if list(positions.values()).count(v) > 1}
            if terminals_pos[n] in duplicate_nodes.values():
                tg.remove_node(n)
                count += 1

    # print(count)
    # print(len(degree_one_nodes(tg)))
    # print(len(tg.nodes))

    if simplify or refine:
        elaborated_lines_dict = {}
        for e in tg.edges:
            name = 'line_' + '_'.join(e)
            elaborated_lines_dict[name] = {}
            elaborated_lines_dict[name]['kwargs'] = {}
            elaborated_lines_dict[name]['Buses'] = [*e]
            elaborated_lines_dict[name]['kwargs']['units'] = tg.edges[e]['units']
            elaborated_lines_dict[name]['kwargs']['length'] = tg.edges[e]['length']
            elaborated_lines_dict[name]['kwargs']['phases'] = tg.edges[e]['phases']
            elaborated_lines_dict[name]['LineCode'] = tg.edges[e]['linecode']
    else:
        elaborated_lines_dict = lines_dict

    # </editor-fold>

    print('Adding lines...')
    for lname, ldata in elaborated_lines_dict.items():
        eulv[tuple(ldata['Buses'])] << kp.Line(**ldata['kwargs']).aka(lname) * lc_el[ldata['LineCode']]

    if insert_loads:
        print('Adding loadprofiles...')
        for lpdata in lp_dict.values():
            eulv << lpdata

    # loads
    td = {'A': '1', 'B': '2', 'C': '3'}
    print('Extracting loads...')
    loads_dict = dict()
    with open(os.path.join(eulv_root, 'Loads.csv'), 'r') as ldf:
        rd = csv.reader(ldf)
        next(rd)  # first line is a comment
        next(rd)  # second line is a comment
        hd = next(rd)  # header without name
        while True:
            try:
                rl = next(rd)
                raw_load_dict = dict(zip(hd, rl))
                loads_dict[raw_load_dict['Name']] = {'Bus': raw_load_dict['Bus'] + '.' + td[raw_load_dict['phases']],
                                                     'duty_name': int(raw_load_dict['Yearly'].split('_')[1]),
                                                     'kwargs': {
                                                     'phases': int(raw_load_dict['numPhases']),
                                                     'conn': raw_load_dict['Connection'],
                                                     'pf': float(raw_load_dict['PF']),
                                                     'kW': float(raw_load_dict['kW']),
                                                     'kV': float(raw_load_dict['kV']),
                                                     'model': raw_load_dict['Model']}}

            except StopIteration:
                break

    if insert_loads:
        print('Adding loads...')
        for ldname, lddata in loads_dict.items():
            if lddata['Bus'].split('.')[0] in eulv.graph().nodes:
                eulv[(lddata['Bus'],)] << kp.Load(**lddata['kwargs'])(model=loadmodel).aka(ldname) * lp_dict[int(lddata['duty_name'])]
            else:
                print('skipped load @refined bus {}'.format(lddata['Bus']))
        noload = ''
    else:
        noload = '_noload'

    print('Linking coordinates...')
    eulv.link_coords(os.path.join(eulv_root, 'Buscoords.csv'))

    this_dir = os.path.dirname(os.path.realpath(__file__))

    if save_ckt:
        eulv.pack_ckt(os.path.join(this_dir, tname + noload + '.zip'))

    if save_json:
        eulv.save_json(os.path.join(this_dir, tname + noload + '.json'))

    if save_effnodes:
        import json
        effnodes = [lddata['Bus'].split('.')[0] for lddata in loads_dict.values() if lddata['Bus'].split('.')[0] in tg.nodes]
        with open(os.path.join(this_dir, 'effective_nodes_' + tname + noload + '.json'), 'w+') as f:
            json.dump(effnodes, f)

    return eulv


if __name__ == '__main__':
    eulv = make_circuit('eu_lv_refined_ns',
                        save_json=True,
                        save_ckt=True,
                        save_effnodes=False,
                        insert_loads=False,
                        simplify=False,
                        refine=True
                        )
    a = 3
    pass
