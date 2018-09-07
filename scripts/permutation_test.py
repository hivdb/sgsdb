#! /usr/bin/env python

import os
import re
import sys
import csv
from hivdbql import app

from collections import Counter

db = app.db
models = app.models

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
DB_AA_VARIANTS_TABLE = os.path.join(
    BASEDIR, 'data', 'prevalence', 'dbAminoAcidVariantsAll.csv')
MUTPATTERN = re.compile(r'^([A-Z])(\d+)([A-Z_*-]+)')


def parse_mutation(mut):
    match = MUTPATTERN.match(mut)
    try:
        cons, pos, aas = match.groups()
    except AttributeError:
        return '', -1, ''
    if '_' in aas:
        # we only care about if it's an insertion
        aas = '_'
    return cons, int(pos), aas


def unusual_mutation_map():
    uum = set()
    with open(DB_AA_VARIANTS_TABLE) as fp:
        data = csv.DictReader(fp)
        for row in data:
            if row['isUsual'] == 'True':
                continue
            uum.add((row['gene'], int(row['position']), row['aa']))
    return uum


UUM = unusual_mutation_map()


def get_random_patients(size, gene):
    Patient = models.Patient
    Isolate = models.Isolate
    ClinicalIsolate = models.ClinicalIsolate
    ptids = (
        db.session.query(Patient.id)
        .filter(
            Patient.id != 11630,
            Patient.isolates.any(db.and_(
                Isolate.gene == gene,
                Isolate.isolate_type == 'Clinical',
                Isolate.clinical_isolate.has(
                    ClinicalIsolate.source == 'Plasma'
                )
            ))
        )
        .order_by(db.func.random())
        .limit(size)
        .all())
    ptids = [i for i, in ptids]
    return Patient.query.filter(Patient.id.in_(ptids))


def count_isolates(patients, gene):
    return sum(len([i for i in pt.isolates if i.gene == gene])
               for pt in patients)


def count_unusual_mutations(gene, muts):
    uum = set()
    for pos, aa in muts:
        if (gene, pos, aa) in UUM:
            uum.add((gene, pos, aa))
    return len(uum)


def parse_mutations(muts):
    result = []
    for mut in muts:
        cons, pos, aas = parse_mutation(mut)
        if pos == -1:
            continue
        for aa in aas:
            if aa == cons:
                continue
            result.append((pos, aa))
    return result


def count_mutations(patients, gene, should_remove_mixtures):
    total = Counter()
    for pt in patients:
        for iso in pt.isolates:
            if iso.gene != gene:
                continue
            for seq in iso.sequences:
                muts = seq.sierra_mutations
                muts = parse_mutations(muts)
                if gene == 'RT':
                    muts = [(p, a) for p, a in muts if p <= 240]
                for mut in muts:
                    total[mut] += 1
    num_muts = len(total.keys())
    num_uums = count_unusual_mutations(gene, total.keys())
    muts_o1 = [k for k, v in total.items() if v == 1]
    num_muts_ge2 = num_muts - len(muts_o1)
    num_uums_ge2 = num_uums - count_unusual_mutations(gene, muts_o1)

    return (
        num_muts, num_uums, num_uums / num_muts,
        num_muts_ge2, num_uums_ge2, num_uums_ge2 / num_muts_ge2
    )


def main():
    size, gene, times = sys.argv[1:]
    print('# Patients', '# Isolates',
          '# Mutations GE1', '# Unusual Mutations GE1',
          '% Unusual Mutations GE1',
          '# Mutations GE2', '# Unusual Mutations GE2',
          '% Unusual Mutations GE2',
          sep='\t')
    for _ in range(int(times)):
        patients = get_random_patients(size, gene).all()
        row = [len(patients)]
        row.append(count_isolates(patients, gene))
        row.extend(count_mutations(patients, gene, False))
        print(*row, sep='\t')
        sys.stdout.flush()


if __name__ == '__main__':
    main()
