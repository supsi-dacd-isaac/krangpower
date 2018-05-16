# this is a function meant to be run una-tantum in order to get all the interface identifiers from the documentation of
# opendssdirect.py and it is not meant to be run dynamically.


def get_methods_from_docs(doc_dir):
    
    import re
    import os
    
    totdic = {}

    for filename in os.listdir(doc_dir):
        if filename.endswith('rst'):
            with open(os.path.join(doc_dir, filename), 'r', encoding='utf8') as f:
                try:
                    r = f.read()
                except UnicodeDecodeError as e:
                    print(filename)
                    print(e)
                    continue
                t = re.findall('(?<=\.\.\ function\:\:\ )(?:[a-z]|[A-Z]|[0-9]|\.)+(?=\n)', r)

                for fcn in t:
                    parts = fcn.split('.')
                    try:
                        totdic[parts[1]].append(parts[2])
                    except KeyError:
                        totdic[parts[1]] = [parts[2]]

    # str(totdic).replace('], ', '],\n').replace("':", '.stack):').replace("\n'", '\ntuple(p_odr.')
    
    return totdic
