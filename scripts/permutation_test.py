#! /usr/bin/env python

import os
import re
import sys
import json
from hivdbql import app

from collections import Counter

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


def get_random_patients(size, gene, profile):
    Patient = models.Patient
    Isolate = models.Isolate
    Sequence = models.Sequence
    ClinicalIsolate = models.ClinicalIsolate
    query = (
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
    )
    ptids = []
    partialsizes = []
    for cat1, criterion1 in SUBTYPE_CATEGORIES.items():
        ratio = profile['{}Ratio'.format(cat1)]
        catsize1 = size * ratio
        for cat2, criterion2 in RX_CATEGORIES.items():
            ratio = profile['{}Ratio'.format(cat2)]
            catsize2 = int(catsize1 * ratio)
            pptids = (
                query.filter(criterion1, criterion2)
                .order_by(db.func.random())
                .limit(catsize2).all())
            partialsizes.append(len(pptids))
            ptids += pptids
    ptids = [i for i, in ptids]
    return (
        Patient.query.filter(Patient.id.in_(ptids))
        .options(
            db.selectinload(Patient.isolates)
            .selectinload(Isolate.sequences)
            .joinedload(Sequence.derived_mutations)),
        partialsizes)


def get_single_isolates(patients, gene):
    isolates = []
    for patient in patients:
        for iso in patient.isolates:
            if iso.gene == gene:
                isolates.append(iso)
                break
    return isolates


def count_unusual_mutations(gene, muts):
    uum = set()
    for pos, aa in muts:
        if (gene, pos, aa) in UUM:
            uum.add((gene, pos, aa))
    return len(uum)


def count_apobec_mutations(gene, muts):
    apm = set()
    for pos, aa in muts:
        if (gene, pos, aa) in APM:
            apm.add((gene, pos, aa))
    return len(apm)


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
    total = Counter()
    for iso in isolates:
        for seq in iso.sequences:
            muts = seq.sierra_mutations
            muts = parse_mutations(muts)
            if iso.gene == 'RT':
                muts = [(p, a) for p, a in muts if p <= 240]
            for mut in muts:
                total[mut] += 1
    num_muts = len(total.keys())
    num_uums = count_unusual_mutations(iso.gene, total.keys())
    num_apms = count_apobec_mutations(iso.gene, total.keys())

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
        patients = patients.all()
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
