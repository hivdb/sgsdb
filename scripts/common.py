import os
import csv
import json
import requests
from decimal import Decimal
from functools import cache  # require Python 3.9
from collections import Counter

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)

FACTSHEET = os.path.join(BASEDIR, 'data', 'SGS.sequences.fact.csv')
SIEERAREPORT = os.path.join(BASEDIR, 'local', 'SGS.sequences.json')
DB_AA_VARIANTS_TABLE = (
    'https://raw.githubusercontent.com/hivdb/hivfacts/'
    'master/data/aapcnt/rx-all_subtype-all.json'
)
APOBEC_TABLE = (
    'https://raw.githubusercontent.com/hivdb/hivfacts/'
    'master/data/apobecs/apobecs.json'
)
AGG_MUTATIONS = {
    'PR': os.path.join(BASEDIR, 'data', 'prevalence', 'CompPR{}.csv'),
    'RT': os.path.join(BASEDIR, 'data', 'prevalence', 'CompRT{}.csv'),
    'IN': os.path.join(BASEDIR, 'data', 'prevalence', 'CompIN{}.csv'),
}
PREC3 = Decimal('1.000')


@cache
def unusual_mutation_map():
    uum = set()
    resp = requests.get(DB_AA_VARIANTS_TABLE)
    data = resp.json()
    for row in data:
        if not row['isUnusual']:
            continue
        uum.add((row['gene'], row['position'], row['aa']))
    return uum


@cache
def apobec_mutation_map():
    resp = requests.get(APOBEC_TABLE)
    return {(r['gene'], r['position'], r['aa']) for r in resp.json()}


def load_sequences(filtered=False):
    sierra_reports = load_sierra_reports()
    with open(FACTSHEET) as fp:
        if fp.read(1) != '\ufeff':
            fp.seek(0)
        result = []
        for seq in csv.DictReader(fp):
            if seq['_Include'] == 'TRUE':
                seq['_Sierra'] = sierra_reports[seq['Accession']]
                result.append(seq)
        counter = Counter()
        for seq in result:
            counter[(seq['PtIdentifier'], seq['CollectionDate'])] += 1
        for seq in result:
            seq['_Filtered'] = \
                counter[(seq['PtIdentifier'], seq['CollectionDate'])] >= 10
        if filtered:
            result = [seq for seq in result if seq['_Filtered']]
        return result


def load_sierra_reports():
    with open(SIEERAREPORT) as fp:
        sequences = json.load(fp)
        return {s['inputSequence']['header'].split('.', 1)[0]: s
                for s in sequences}


def load_aggregated_mutations(gene, subset='All'):
    with open(AGG_MUTATIONS[gene].format(subset)) as fp:
        result = []
        for row in csv.DictReader(fp):
            row['excluded'] = False
            if row['AA'] in ('-', '_', '*', 'X'):
                row['excluded'] = True
            if row['dbPcnt'] == 'NA':
                row['dbPcnt'] = .0
                row['pcntFold'] = 0xffff
            for k in ('Pos', 'Count', 'PosTotal',
                      'PatientCount', 'PatientPosTotal',
                      'SampleCount', 'SamplePosTotal'):
                row[k] = int(row[k])
            for k in ('sgsPcnt', 'dbPcnt', 'pcntFold'):
                row[k] = float(row[k])
            for k in ('IsAPOBEC', 'isUnusual'):
                row[k] = row[k].upper() == 'TRUE'
            result.append(row)
        return result
