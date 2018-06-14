import csv
import os
import requests
import io
import zipfile
import krangpower as kp

um = kp.UM

test_dir = os.path.join(os.getenv("AppData"), 'krangpower//test')


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
            lines_dict[rl[0]]['kwargs']['length'] = float(ll['Length'])
            lines_dict[rl[0]]['kwargs']['phases'] = len(ll['Phases'])
            lines_dict[rl[0]]['LineCode'] = ll['LineCode'].replace('.', '')
        except StopIteration:
            break

td = {'A': '1', 'B': '2', 'C': '3'}

# loads
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

# the data is contained in:
# skw, loads_dict, lines_dict, lc_dict, lp_dict, transformer
print('Creating circuit...')
SRC = kp.Vsource(**skw)
eulv = kp.Krang('eu_lv_test', SRC)

eulv.set(basefreq=50.0 * um.Hz, voltagebases=[11.0, 0.416] * um.kV)
eulv.command('calcvoltagebases')

print('Adding transformer...')
eulv['sourcebus', '1'] << TR1.aka('TR1')

print('Adding linecodes...')
lc_el = {}
for lcname, lcdata in lc_dict.items():
    lc_el[lcname] = kp.LineCode_S(lcname, **lcdata)
    eulv << lc_el[lcname]

print('Adding loadprofiles...')
for lpdata in lp_dict.values():
    eulv << lpdata

print('Adding lines...')
for lname, ldata in lines_dict.items():
    eulv[tuple(ldata['Buses'])] << kp.Line(**ldata['kwargs']).aka(lname) * lc_el[ldata['LineCode']]

print('Adding loads...')
for ldname, lddata in loads_dict.items():
    eulv[(lddata['Bus'],)] << kp.Load(**lddata['kwargs']).aka(ldname) * lp_dict[int(lddata['duty_name'])]

eulv.link_coords(os.path.join(eulv_root, 'Buscoords.csv'))

eulv.set(basefreq=50.0 * um.Hz)

eulv.set(time=[0, 0], number=1, stepsize=1.0 * um.min)
eulv.solve()
p_1 = eulv.export('powers')
print(p_1)

eulv.set(time=[0, 60*565], number=1, stepsize=1.0 * um.min)
eulv.solve()
p_566 = eulv.export('powers')
print(p_566)

eulv.set(time=[0, 60*1439], number=1, stepsize=1.0 * um.min)
eulv.solve()
v_1440 = eulv.export('voltages')
print(v_1440)

eulv.set(time=[0, 0], number=1440, stepsize=1.0 * um.min)
eulv << kp.Monitor().aka('m1') * eulv['load.LOAD1']
eulv << kp.Monitor().aka('m23') * eulv['load.LOAD32']
eulv << kp.Monitor().aka('m53') * eulv['load.LOAD53']
eulv << kp.Monitor().aka('mtr') * eulv['transformer.TR1']
eulv.solve()

m1_out = eulv.export('monitor m1')
m23_out = eulv.export('monitor m23')
m53_out = eulv.export('monitor m53')
mtr_out = eulv.export('monitor mtr')

print(m1_out)
print(m23_out)
print(m53_out)
print(mtr_out)

eulv.pack_ckt(r'D:/GDrive/Pycharm/krangpower/test/eu_lv/eu_lv.zip')
eulv.save_json(r'D:/GDrive/Pycharm/krangpower/test/eu_lv/eu_lv.json')
