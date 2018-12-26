#! /usr/bin/env python
import os
import csv

from decimal import Decimal
from collections import Counter, defaultdict

from numpy import percentile
from scipy.stats import linregress

from common import (load_sequences, load_aggregated_mutations,
                    apobec_mutation_map, BASEDIR, PREC3)

REPORT_PATH = os.path.join(BASEDIR, 'data', 'report.csv')
GENES = ('PR', 'RT', 'IN')
TREATMENTS = ('ART', 'None', 'Unknown')
CATEGORIES = ('All', 'SubtypeB', 'SubtypeC', 'Non-SubtypeBC', 'Naive', 'ART')
GENE_RANGES = {
    'PR': [
        ('FirstAA=1, LastAA=99', lambda f, l: f == 1 and l == 99),
    ],
    'RT': [
        ('FirstAA=1, LastAA>=240', lambda f, l: f == 1 and l >= 240),
        ('FirstAA=1, LastAA=560', lambda f, l: f == 1 and l >= 560),
    ],
    'IN': [
        ('FirstAA=1, LastAA=288', lambda f, l: f == 1 and l == 288),
    ],
}
SUBTYPES = ('B', 'C', 'Other')
APM = apobec_mutation_map()

header = ['name', 'subset', 'value', 'percent',
          'range_0', 'range_100', 'r_squared',
          'p_value', 'stderr']


def make_row(name, subset, value, **kws):
    kws.update({
        'name': name,
        'subset': subset,
        'value': value
    })
    if 'total' in kws:
        kws['percent'] = '{}%'.format(Decimal(
            kws['value'] / kws['total'] * 100
        ).quantize(PREC3))
    row = {}
    for h in header:
        row[h] = kws.get(h)
    return row


def make_percentile_row(name, subset, arr):
    arr = list(arr)
    p0, median, p100 = percentile(arr, (0, 50, 100))
    kws = {
        'name': name,
        'subset': subset,
        'value': median,
        'range_0': p0,
        'range_100': p100,
    }
    return make_row(**kws)


def make_linregress_row(name, subset, arr1, arr2=None):
    arr1 = list(arr1)
    arr2 = None if arr2 is None else list(arr2)
    _, _, r_value, p_value, stderr = linregress(arr1, arr2)
    kws = {
        'name': name,
        'subset': subset,
        'value': None,
        'r_squared': r_value ** 2,
        'p_value': p_value,
        'stderr': stderr
    }
    return make_row(**kws)


def get_subtype(seq):
    subtype = seq['subtypeText'].split(' (', 1)[0]
    if subtype in SUBTYPES:
        return subtype
    else:
        return 'Other'


def prevalence_stat(gene, cat, ptcount):
    mutations = load_aggregated_mutations(gene, cat)
    ptmuttotal = sum(m['PatientCount'] for m in mutations
                     if not m['excluded']) / ptcount
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}'.format(gene, cat),
        ptmuttotal)
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}, '
        'IsUnusual, NotAPOBEC'.format(gene, cat),
        sum(m['PatientCount'] for m in mutations if
            not m['isUsual'] and not m['IsAPOBEC'] and
            not m['excluded']) / ptcount,
        total=ptmuttotal)
    ptseqtotal = sum(m['PatientCount'] for m in mutations
                     if m['Count'] > 1 and not m['excluded']) / ptcount
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}, NumSequences>1'.format(gene, cat),
        ptseqtotal, total=ptmuttotal)
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}, NumSequences>1, '
        'IsUnusual, NotAPOBEC'.format(gene, cat),
        sum(m['PatientCount'] for m in mutations if m['Count'] > 1
            and not m['isUsual'] and not m['IsAPOBEC'] and
            not m['excluded']) / ptcount,
        total=ptseqtotal)
    ptpttotal = sum(m['PatientCount'] for m in mutations
                    if m['PatientCount'] > 1 and not m['excluded']) / ptcount
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}, NumPatients>1'.format(gene, cat),
        ptpttotal, total=ptmuttotal)
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}, NumPatients>1, '
        'IsUnusual, NotAPOBEC'.format(gene, cat),
        sum(m['PatientCount'] for m in mutations if m['PatientCount'] > 1
            and not m['isUsual'] and not m['IsAPOBEC'] and
            not m['excluded']) / ptcount,
        total=ptpttotal)
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}, NumSequences>1, '
        'IsAPOBEC'.format(gene, cat),
        sum(m['PatientCount'] for m in mutations if m['Count'] > 1
            and m['IsAPOBEC']) / ptcount)
    yield make_row(
        '# Mutations per Patient',
        'Gene={}, Category={}, NumPatients>1, '
        'IsAPOBEC'.format(gene, cat),
        sum(m['PatientCount'] for m in mutations if m['PatientCount'] > 1
            and m['IsAPOBEC']) / ptcount)

    muttotal = len([m for m in mutations if not m['excluded']])
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}'.format(gene, cat),
        muttotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, '
        'IsUnusual, NotAPOBEC'.format(gene, cat),
        len([m for m in mutations if
             not m['isUsual'] and not m['IsAPOBEC'] and not m['excluded']]),
        total=muttotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, '
        'IsAPOBEC'.format(gene, cat),
        len([m for m in mutations if m['IsAPOBEC']]))
    seqtotal = len([m for m in mutations
                    if m['Count'] > 1 and not m['excluded']])
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, NumSequences>1'.format(gene, cat),
        seqtotal, total=muttotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, NumSequences>1, '
        'IsUnusual, NotAPOBEC'.format(gene, cat),
        len([m for m in mutations if m['Count'] > 1
            and not m['isUsual'] and not m['IsAPOBEC'] and not m['excluded']]),
        total=seqtotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, NumSequences>1, '
        'IsAPOBEC'.format(gene, cat),
        len([m for m in mutations if m['Count'] > 1 and m['IsAPOBEC']]))
    pttotal = len([m for m in mutations
                   if m['PatientCount'] > 1 and not m['excluded']])
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, NumPatients>1'.format(gene, cat),
        pttotal, total=muttotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, NumPatients>1, '
        'IsUnusual, NotAPOBEC'.format(gene, cat),
        len([m for m in mutations if m['PatientCount'] > 1
             and not m['isUsual'] and not m['IsAPOBEC']
             and not m['excluded']]),
        total=pttotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, NumPatients>1, '
        'IsAPOBEC'.format(gene, cat),
        len([m for m in mutations if m['PatientCount'] > 1 and m['IsAPOBEC']]))
    yield make_linregress_row(
        'Prevalence Correlation b/t SGS and HIVDB',
        'Gene={}, Category={}'.format(gene, cat),
        [(m['sgsPcnt'], m['dbPcnt']) for m in mutations if not m['excluded']])


def overall_prevalence_stat():
    mutations = sum([
        load_aggregated_mutations(gene, 'All') for gene in GENES
    ], [])
    muttotal = len(mutations)
    pttotal = len([m for m in mutations if m['PatientCount'] > 1])
    yield make_row(
        '# Uniq. Mutations',
        'Category=All, NumPatients>1',
        pttotal, total=muttotal)
    pttotal2 = len([m for m in mutations if m['PatientCount'] > 2])
    yield make_row(
        '# Uniq. Mutations',
        'Category=All, NumPatients>2',
        pttotal2, total=muttotal)


def basic_stat(sequences):
    pmids = set()
    ptids = set()
    pttps = defaultdict(set)
    geneseqs = Counter()
    generangeseqs = Counter()
    geneptids = defaultdict(set)
    genepttps = defaultdict(set)
    rxptids = defaultdict(set)
    subtypeseqs = Counter()
    apobecseqs = Counter()
    subtypeptids = defaultdict(set)
    pttpseqs = Counter()
    for seq in sequences:
        pmids.add(seq['MedlineID'])
        sierra = seq['_Sierra']
        ptid = seq['PtIdentifier']
        pttp = (ptid, seq['CollectionDate'])
        ptids.add(ptid)
        rx = seq['Rx']
        rxptids[rx].add(ptid)
        pttps[ptid].add(pttp)
        pttpseqs[pttp] += 1
        subtype = get_subtype(sierra)
        subtypeseqs[subtype] += 1
        subtypeptids[subtype].add(ptid)
        for geneseq in sierra['alignedGeneSequences']:
            gene = geneseq['gene']['name']
            geneseqs[gene] += 1
            geneptids[gene].add(ptid)
            subtypeptids[(subtype, gene)].add(ptid)
            rxptids[(rx, gene)].add(ptid)
            genepttps[gene].add(pttp)
            for (subset, func) in GENE_RANGES[gene]:
                if func(geneseq['firstAA'], geneseq['lastAA']):
                    generangeseqs['Gene={}, {}'.format(gene, subset)] += 1
            numapobecs = 0
            for mut in geneseq['mutations']:
                pos = mut['position']
                if gene == 'RT' and pos > 240:
                    continue
                aa = mut['AAs']
                if (gene, pos, aa) in APM:
                    numapobecs += 1
            numapobecs = min(3, numapobecs)
            if numapobecs > 0:
                apobecseqs[(gene, numapobecs)] += 1

    yield make_row('# Studies', None, len(pmids))

    totalseqs = len(sequences)
    yield make_row('# Sequences', None, totalseqs)
    for gene in GENES:
        yield make_row('# Sequences',
                       'Gene={}'.format(gene),
                       geneseqs[gene],
                       total=totalseqs)
        addcond = ', MinPos=1, MaxPos=240' if gene == 'RT' else ''
        yield make_row('# Sequences',
                       'Gene={}{}, NumAPOBECs=1'.format(gene, addcond),
                       apobecseqs[(gene, 1)],
                       total=totalseqs)
        yield make_row('# Sequences',
                       'Gene={}{}, NumAPOBECs=2'.format(gene, addcond),
                       apobecseqs[(gene, 2)],
                       total=totalseqs)
        yield make_row('# Sequences',
                       'Gene={}{}, NumAPOBECs>=3'.format(gene, addcond),
                       apobecseqs[(gene, 3)],
                       total=totalseqs)

    for subset, value in generangeseqs.items():
        yield make_row('# Sequences', subset, value, total=totalseqs)

    for subtype in SUBTYPES:
        yield make_row('# Sequences',
                       'Subtype={}'.format(subtype),
                       subtypeseqs[subtype],
                       total=totalseqs)

    totalpts = len(ptids)
    yield make_row('# Patients', None, totalpts)
    for gene in GENES:
        yield make_row('# Patients',
                       'Gene={}'.format(gene),
                       len(geneptids[gene]),
                       total=totalpts)
    for subtype in SUBTYPES:
        yield make_row('# Patients',
                       'Subtype={}'.format(subtype),
                       len(subtypeptids[subtype]),
                       total=totalpts)
        for gene in GENES:
            yield make_row('# Patients',
                           'Subtype={}, Gene={}'.format(subtype, gene),
                           len(subtypeptids[(subtype, gene)]),
                           total=totalpts)
    for rx in TREATMENTS:
        yield make_row('# Patients',
                       'Rx={}'.format(rx),
                       len(rxptids[rx]),
                       total=totalpts)
        for gene in GENES:
            yield make_row('# Patients',
                           'Rx={}, Gene={}'.format(rx, gene),
                           len(rxptids[(rx, gene)]),
                           total=totalpts)
    yield make_row('# Patients', 'NumSample=1',
                   len([k for k, v in pttps.items() if len(v) == 1]),
                   total=totalpts)
    yield make_row('# Patients', 'NumSample>1',
                   len([k for k, v in pttps.items() if len(v) > 1]),
                   total=totalpts)

    totalpttps = sum(len(v) for v in pttps.values())
    yield make_row('# Samples (Patient Time Points)', None, totalpttps)
    for gene in GENES:
        yield make_row('# Samples (Patient Time Points)',
                       'Gene={}'.format(gene),
                       len(genepttps[gene]),
                       total=totalpttps)

    yield make_percentile_row("Med. Samples per Patient",
                              None, [len(v) for v in pttps.values()])
    yield make_percentile_row("Med. Sequences per Sample",
                              None, pttpseqs.values())


def main():
    sequences = load_sequences()
    with open(REPORT_PATH, 'w') as fp:
        writer = csv.DictWriter(fp, header)
        writer.writeheader()
        ptcount = {}
        for row in basic_stat(sequences):
            writer.writerow(row)
            if row['name'] == '# Patients':
                if row['subset'] == 'Gene=PR':
                    ptcount['PR'] = row['value']
                elif row['subset'] == 'Gene=RT':
                    ptcount['RT'] = row['value']
                elif row['subset'] == 'Gene=IN':
                    ptcount['IN'] = row['value']
        for cat in CATEGORIES:
            for gene in GENES:
                writer.writerows(prevalence_stat(gene, cat, ptcount[gene]))
        writer.writerows(overall_prevalence_stat())


if __name__ == '__main__':
    main()
