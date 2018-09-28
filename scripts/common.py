import os
import csv
import json
from decimal import Decimal

from hivdbql import app

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)

FACTSHEET = os.path.join(BASEDIR, 'data', 'SGS.sequences.fact.csv')
SIEERAREPORT = os.path.join(BASEDIR, 'local', 'SGS.sequences.json')
DB_AA_VARIANTS_TABLE = os.path.join(
    BASEDIR, 'data', 'prevalence', 'dbAminoAcidVariantsAll.csv')

AGG_MUTATIONS = {
    'PR': os.path.join(BASEDIR, 'data', 'prevalence', 'CompPR{}.csv'),
    'RT': os.path.join(BASEDIR, 'data', 'prevalence', 'CompRT{}.csv'),
    'IN': os.path.join(BASEDIR, 'data', 'prevalence', 'CompIN{}.csv'),
}
PREC3 = Decimal('1.000')

models = app.models


def unusual_mutation_map():
    uum = set()
    with open(DB_AA_VARIANTS_TABLE) as fp:
        data = csv.DictReader(fp)
        for row in data:
            if row['isUsual'] == 'True':
                continue
            uum.add((row['gene'], int(row['position']), row['aa']))
    return uum


def apobec_mutation_map():
    apobecs = models.LUAPOBEC.query.filter(
        models.LUAPOBEC.is_apobec.is_(True)
    )
    return {(m.gene, m.pos, m.hm) for m in apobecs}


def load_sequences():
    sierra_reports = load_sierra_reports()
    with open(FACTSHEET) as fp:
        if fp.read(1) != '\ufeff':
            fp.seek(0)
        result = []
        for seq in csv.DictReader(fp):
            if seq['_Include'] == 'TRUE':
                seq['_Sierra'] = sierra_reports[seq['Accession']]
                result.append(seq)
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
            for k in ('Pos', 'Count', 'PosTotal',
                      'PatientCount', 'PatientPosTotal'):
                row[k] = int(row[k])
            for k in ('sgsPcnt', 'dbPcnt', 'pcntFold'):
                row[k] = float(row[k])
            for k in ('isAPOBEC', 'isUsual'):
                row[k] = row[k] == 'True'
            result.append(row)
        return result
