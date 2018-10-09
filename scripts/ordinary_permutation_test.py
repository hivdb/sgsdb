#! /usr/bin/env python

import os
import re
import sys
import json
import random
from hivdbql import app

from common import unusual_mutation_map, apobec_mutation_map

db = app.db
models = app.models

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
PROFILE_PATH = os.path.join(BASEDIR, 'local', 'permutation_profile.json')
MUTPATTERN = re.compile(r'^([A-Z])(\d+)([A-Z_*-]+)')

UUM = unusual_mutation_map()
APM = apobec_mutation_map()

SUBTYPE_CATEGORIES = {
    'SubtypeB': models.Patient.isolates.any(
        models.Isolate._subtype.has(
            models.Subtype.subtype == 'B'
        ),
    ),
    'SubtypeC': models.Patient.isolates.any(
        models.Isolate._subtype.has(
            models.Subtype.subtype == 'C'
        ),
    ),
    'SubtypeOther': ~models.Patient.isolates.any(
        models.Isolate._subtype.has(
            models.Subtype.subtype.in_(['B', 'C'])
        ),
    )
}
RX_CATEGORIES = {
    'RxART': models.Patient.treatments.any(
        models.RxHistory.regimen_name != 'None'
    ),
    'RxNaive': models.Patient.treatments.any(
        models.RxHistory.regimen_name == 'None'
    ),
}

PATIENT_POOL = {
    (c1, c2): []
    for c1 in SUBTYPE_CATEGORIES
    for c2 in RX_CATEGORIES
}


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


def load_all(query, limit=5000):
    offset = 0
    while True:
        result = query.offset(offset).limit(limit).all()
        yield from result
        offset += limit
        print(offset, len(result), file=sys.stderr)
        if len(result) < limit:
            break


def patient_pool(cat1, cat2, gene):
    if not PATIENT_POOL[(cat1, cat2)]:
        Patient = models.Patient
        Isolate = models.Isolate
        Sequence = models.Sequence
        ClinicalIsolate = models.ClinicalIsolate
        query = (
            Patient.query
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
            .options(
                db.selectinload(Patient.isolates)
                .selectinload(Isolate.sequences)
                .joinedload(Sequence.derived_mutations))
        )
        for c1, criterion1 in SUBTYPE_CATEGORIES.items():
            for c2, criterion2 in RX_CATEGORIES.items():
                print(c1, c2, file=sys.stderr)
                parquery = (
                    query.filter(criterion1, criterion2)
                    .order_by(Patient.id))
                PATIENT_POOL[(c1, c2)] = list(load_all(parquery))
    return PATIENT_POOL[(cat1, cat2)]


def get_random_patients(size, gene, profile):
    patients = []
    partialsizes = []
    for cat1, criterion1 in SUBTYPE_CATEGORIES.items():
        ratio = profile['{}Ratio'.format(cat1)]
        catsize1 = size * ratio
        for cat2, criterion2 in RX_CATEGORIES.items():
            ratio = profile['{}Ratio'.format(cat2)]
            catsize2 = int(catsize1 * ratio)
            pool = patient_pool(cat1, cat2, gene)
            ppatients = random.sample(pool, catsize2)
            partialsizes.append(len(ppatients))
            patients += ppatients
    return patients, partialsizes


def get_single_isolates(patients, gene):
    isolates = []
    for patient in patients:
        for iso in patient.isolates:
            if iso.gene == gene:
                isolates.append(iso)
                break
    return isolates


def count_unusual_mutations(gene, muts):
    uum = 0
    for pos, aa in muts:
        if (gene, pos, aa) in UUM:
            uum += 1
    return uum


def count_apobec_mutations(gene, muts):
    apm = 0
    for pos, aa in muts:
        if (gene, pos, aa) in APM:
            apm += 1
    return apm


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


def count_mutations(isolates):
    total = []
    for iso in isolates:
        for seq in iso.sequences:
            muts = seq.sierra_mutations
            muts = parse_mutations(muts)
            if iso.gene == 'RT':
                muts = [(p, a) for p, a in muts if p <= 240]
            for mut in muts:
                total.append(mut)
    num_muts = len(total)
    num_uums = count_unusual_mutations(iso.gene, total)
    num_apms = count_apobec_mutations(iso.gene, total)

    return (
        num_muts, num_uums, num_uums / num_muts,
        num_apms, num_apms / num_muts
    )


def entrypoint(gene, profile, times):
    size = profile['{}NumPatients'.format(gene)]
    print('# Patients', '# Isolates',
          '# Mutations', '# Unusual Mutations',
          '% Unusual Mutations',
          '# APOBEC Mutations',
          '% APOBEC Mutations',
          *['# Pts ({} {})'.format(s, r)
            for s in SUBTYPE_CATEGORIES
            for r in RX_CATEGORIES],
          sep='\t')
    for _ in range(int(times)):
        patients, psizes = get_random_patients(size, gene, profile)
        row = [len(patients)]
        isolates = get_single_isolates(patients, gene)
        row.append(len(isolates))
        row.extend(count_mutations(isolates))
        row.extend(psizes)
        print(*row, sep='\t')
        sys.stdout.flush()


def main():
    if len(sys.argv) != 3:
        print('Usage: {} <GENE> <REPEAT>'
              .format(sys.argv[0]), file=sys.stderr)
        exit(1)
    gene, times = sys.argv[1:]
    with open(PROFILE_PATH) as fp:
        profile = json.load(fp)
    entrypoint(gene, profile, times)


if __name__ == '__main__':
    main()
