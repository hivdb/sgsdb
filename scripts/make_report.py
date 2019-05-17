#! /usr/bin/env python
import os
import csv

from decimal import Decimal
from collections import Counter, defaultdict

from numpy import percentile
from scipy.stats import linregress, chi2_contingency

from common import (load_sequences, load_aggregated_mutations,
                    apobec_mutation_map, BASEDIR, PREC3)

REPORT_PATH = os.path.join(BASEDIR, 'data', 'report.csv')
GENES = ('PR', 'RT', 'IN')
TREATMENTS = ('ART', 'None', 'Unknown')
CATEGORIES = ('All', 'SubtypeB', 'SubtypeC', 'Non-SubtypeBC', 'Naive', 'ART')
GENE_RANGES = {
    'PR': [
        ('FirstAA=1, LastAA=99', lambda fa, la: fa == 1 and la == 99),
    ],
    'RT': [
        ('FirstAA=1, LastAA>=240', lambda fa, la: fa == 1 and la >= 240),
        ('FirstAA=1, LastAA=560', lambda fa, la: fa == 1 and la >= 560),
    ],
    'IN': [
        ('FirstAA=1, LastAA=288', lambda fa, la: fa == 1 and la == 288),
    ],
}
SUBTYPES = ('B', 'C', 'Other')
APM = apobec_mutation_map()

header = ['name', 'subset', 'value', 'percent',
          'range_0', 'range_100', 'r_squared',
          'p_value', 'stderr', 'note']


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
    p0, p25, median, p75, p100 = percentile(arr, (0, 25, 50, 75, 100))
    mean = sum(arr) / len(arr)
    kws = {
        'name': name,
        'subset': subset,
        'value': median,
        'range_0': p0,
        'range_100': p100,
        'note': 'IQR: {} - {}; Mean: {}'.format(p25, p75, mean)
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


def prevalence_stat(gene, cat, splcount):
    mutations = load_aggregated_mutations(gene, cat)
    splmuttotal = sum(m['SampleCount'] for m in mutations
                      if not m['excluded'])
    yield make_row(
        '# Mutations per Sample',
        'Gene={}, Category={}'.format(gene, cat),
        splmuttotal / splcount)
    spluumtotal = sum(m['SampleCount'] for m in mutations if
                      m['isUnusual'] and not m['excluded'])
    yield make_row(
        '# Mutations per Sample',
        'Gene={}, Category={}, '
        'IsUnusual'.format(gene, cat),
        spluumtotal / splcount, total=splmuttotal / splcount)
    splmutgt1total = sum(m['SampleCount'] for m in mutations
                         if m['Count'] > 1 and not m['excluded'])
    yield make_row(
        '# Mutations per Sample',
        'Gene={}, Category={}, NumSequences>1'.format(gene, cat),
        splmutgt1total / splcount, total=splmuttotal / splcount)

    spluumeq1total = sum(
        m['SampleCount'] for m in mutations if m['Count'] == 1
        and m['isUnusual'] and not m['excluded'])
    spluumgt1total = sum(
        m['SampleCount'] for m in mutations if m['Count'] > 1
        and m['isUnusual'] and not m['excluded'])
    yield make_row(
        '# Mutations',
        'Gene={}, Category={}, NumSequences=1, '
        'IsUnusual'.format(gene, cat),
        spluumeq1total,
        total=spluumtotal)
    _, p, _, _ = chi2_contingency([
        [spluumeq1total, splmuttotal - spluumeq1total],
        [spluumgt1total, splmuttotal - spluumgt1total]
    ])
    yield make_row(
        '# Mutations',
        'Gene={}, Category={}, NumSequences>1, '
        'IsUnusual'.format(gene, cat),
        spluumgt1total,
        total=spluumtotal,
        p_value=p)
    yield make_row(
        '# Mutations per Sample',
        'Gene={}, Category={}, NumSequences=1, '
        'IsUnusual'.format(gene, cat),
        spluumeq1total / splcount,
        total=splmuttotal / splcount)
    yield make_row(
        '# Mutations per Sample',
        'Gene={}, Category={}, NumSequences>1, '
        'IsUnusual'.format(gene, cat),
        spluumgt1total / splcount,
        total=splmuttotal / splcount)
    # splpttotal = sum(m['SampleCount'] for m in mutations
    #                  if m['PatientCount'] > 1
    #                  and not m['excluded']) / splcount
    # yield make_row(
    #     '# Mutations per Sample',
    #     'Gene={}, Category={}, NumPatients>1'.format(gene, cat),
    #     splpttotal, total=splmuttotal)
    # yield make_row(
    #     '# Mutations per Sample',
    #     'Gene={}, Category={}, NumPatients>1, '
    #     'IsUnusual'.format(gene, cat),
    #     sum(m['SampleCount'] for m in mutations if m['PatientCount'] > 1
    #         and m['isUnusual'] and
    #         not m['excluded']) / splcount,
    #     total=splpttotal)
    # yield make_row(
    #     '# Mutations per Sample',
    #     'Gene={}, Category={}, NumSequences>1, '
    #     'IsAPOBEC'.format(gene, cat),
    #     sum(m['SampleCount'] for m in mutations if m['Count'] > 1
    #         and m['IsAPOBEC']) / splcount)
    # yield make_row(
    #     '# Mutations per Sample',
    #     'Gene={}, Category={}, NumPatients>1, '
    #     'IsAPOBEC'.format(gene, cat),
    #     sum(m['SampleCount'] for m in mutations if m['PatientCount'] > 1
    #         and m['IsAPOBEC']) / splcount)

    muttotal = len([m for m in mutations if not m['excluded']])
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}'.format(gene, cat),
        muttotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, '
        'IsUnusual'.format(gene, cat),
        len([m for m in mutations if
             m['isUnusual'] and not m['excluded']]),
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
        'IsUnusual'.format(gene, cat),
        len([m for m in mutations if m['Count'] > 1
            and m['isUnusual'] and not m['excluded']]),
        total=seqtotal)
    yield make_row(
        '# Uniq. Mutations',
        'Gene={}, Category={}, NumSequences>1, '
        'IsAPOBEC'.format(gene, cat),
        len([m for m in mutations if m['Count'] > 1 and m['IsAPOBEC']]))
    # pttotal = len([m for m in mutations
    #                if m['PatientCount'] > 1 and not m['excluded']])
    # yield make_row(
    #     '# Uniq. Mutations',
    #     'Gene={}, Category={}, NumPatients>1'.format(gene, cat),
    #     pttotal, total=muttotal)
    # yield make_row(
    #     '# Uniq. Mutations',
    #     'Gene={}, Category={}, NumPatients>1, '
    #     'IsUnusual'.format(gene, cat),
    #     len([m for m in mutations if m['PatientCount'] > 1
    #          and m['isUnusual'] and not m['excluded']]),
    #     total=pttotal)
    # yield make_row(
    #     '# Uniq. Mutations',
    #     'Gene={}, Category={}, NumPatients>1, '
    #     'IsAPOBEC'.format(gene, cat),
    #     len([m for m in mutations if
    #          m['PatientCount'] > 1 and m['IsAPOBEC']]))
    yield make_linregress_row(
        'Prevalence Correlation b/t SGS and HIVDB',
        'Gene={}, Category={}'.format(gene, cat),
        [(m['sgsPcnt'], m['dbPcnt']) for m in mutations if not m['excluded']])


# def overall_prevalence_stat():
#     mutations = sum([
#         load_aggregated_mutations(gene, 'All') for gene in GENES
#     ], [])
#     muttotal = len(mutations)
#     pttotal = len([m for m in mutations if m['PatientCount'] > 1])
#     yield make_row(
#         '# Uniq. Mutations',
#         'Category=All, NumPatients>1',
#         pttotal, total=muttotal)
#     pttotal2 = len([m for m in mutations if m['PatientCount'] > 2])
#     yield make_row(
#         '# Uniq. Mutations',
#         'Category=All, NumPatients>2',
#         pttotal2, total=muttotal)


def basic_stat(sequences, outcond=''):
    pmids = set()
    ptids = set()
    pttps = defaultdict(set)
    geneseqs = Counter()
    generangeseqs = Counter()
    geneptids = defaultdict(set)
    genepttps = defaultdict(set)
    subtypepttps = defaultdict(set)
    rxptids = defaultdict(set)
    rxpttps = defaultdict(set)
    subtypeseqs = Counter()
    apobecseqs = Counter()
    subtypeptids = defaultdict(set)
    pttpseqs = Counter()
    outcond0 = outcond.strip(' ,')
    for seq in sequences:
        pmids.add(seq['MedlineID'])
        sierra = seq['_Sierra']
        ptid = seq['PtIdentifier']
        pttp = (ptid, seq['CollectionDate'])
        ptids.add(ptid)
        rx = seq['Rx']
        rxptids[rx].add(ptid)
        rxpttps[rx].add(pttp)
        pttps[ptid].add(pttp)
        pttpseqs[pttp] += 1
        subtype = get_subtype(sierra)
        subtypeseqs[subtype] += 1
        subtypeptids[subtype].add(ptid)
        subtypepttps[subtype].add(pttp)
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

    yield make_row('# Studies', outcond0, len(pmids))

    totalseqs = len(sequences)
    yield make_row('# Sequences', outcond0, totalseqs)
    for gene in GENES:
        yield make_row('# Sequences',
                       outcond + 'Gene={}'.format(gene),
                       geneseqs[gene],
                       total=totalseqs)
        addcond = ', MinPos=1, MaxPos=240' if gene == 'RT' else ''
        yield make_row('# Sequences',
                       outcond + 'Gene={}{}, NumAPOBECs=1'
                       .format(gene, addcond),
                       apobecseqs[(gene, 1)],
                       total=totalseqs)
        yield make_row('# Sequences',
                       outcond + 'Gene={}{}, NumAPOBECs=2'
                       .format(gene, addcond),
                       apobecseqs[(gene, 2)],
                       total=totalseqs)
        yield make_row('# Sequences',
                       outcond + 'Gene={}{}, NumAPOBECs>=3'
                       .format(gene, addcond),
                       apobecseqs[(gene, 3)],
                       total=totalseqs)

    for subset, value in generangeseqs.items():
        yield make_row('# Sequences', outcond + subset, value, total=totalseqs)

    for subtype in SUBTYPES:
        yield make_row('# Sequences',
                       outcond + 'Subtype={}'.format(subtype),
                       subtypeseqs[subtype],
                       total=totalseqs)

    totalpts = len(ptids)
    yield make_row('# Patients', outcond0, totalpts)
    for gene in GENES:
        yield make_row('# Patients',
                       outcond + 'Gene={}'.format(gene),
                       len(geneptids[gene]),
                       total=totalpts)
    for subtype in SUBTYPES:
        yield make_row('# Patients',
                       outcond + 'Subtype={}'.format(subtype),
                       len(subtypeptids[subtype]),
                       total=totalpts)
        for gene in GENES:
            yield make_row('# Patients',
                           outcond + 'Subtype={}, Gene={}'
                           .format(subtype, gene),
                           len(subtypeptids[(subtype, gene)]),
                           total=totalpts)
    for rx in TREATMENTS:
        yield make_row('# Patients',
                       outcond + 'Rx={}'.format(rx),
                       len(rxptids[rx]),
                       total=totalpts)
        for gene in GENES:
            yield make_row('# Patients',
                           outcond + 'Rx={}, Gene={}'.format(rx, gene),
                           len(rxptids[(rx, gene)]),
                           total=totalpts)
    yield make_row('# Patients', outcond + 'NumSample=1',
                   len([k for k, v in pttps.items() if len(v) == 1]),
                   total=totalpts)
    yield make_row('# Patients', outcond + 'NumSample>1',
                   len([k for k, v in pttps.items() if len(v) > 1]),
                   total=totalpts)

    totalpttps = sum(len(v) for v in pttps.values())
    yield make_row('# Samples (Patient Time Points)', outcond0, totalpttps)
    for gene in GENES:
        yield make_row('# Samples (Patient Time Points)',
                       outcond + 'Gene={}'.format(gene),
                       len(genepttps[gene]),
                       total=totalpttps)
    for subtype in SUBTYPES:
        yield make_row('# Samples (Patient Time Points)',
                       outcond + 'Subtype={}'.format(subtype),
                       len(subtypepttps[subtype]),
                       total=totalpttps)
    for rx in TREATMENTS:
        yield make_row('# Samples (Patient Time Points)',
                       outcond + 'Rx={}'.format(rx),
                       len(rxpttps[rx]),
                       total=totalpttps)

    yield make_percentile_row("Med. Samples per Patient",
                              outcond0, [len(v) for v in pttps.values()])
    yield make_percentile_row("Med. Sequences per Sample",
                              outcond0, pttpseqs.values())


def main():
    sequences = load_sequences()
    filtered_sequences = [s for s in sequences if s['_Filtered']]
    with open(REPORT_PATH, 'w') as fp:
        writer = csv.DictWriter(fp, header)
        writer.writeheader()
        splcount = {}
        for row in basic_stat(sequences):
            writer.writerow(row)

        for row in basic_stat(filtered_sequences, 'SequencesPerSample>9, '):
            writer.writerow(row)
            if row['name'] == '# Samples (Patient Time Points)':
                if row['subset'].endswith('Gene=PR'):
                    splcount['PR'] = row['value']
                elif row['subset'].endswith('Gene=RT'):
                    splcount['RT'] = row['value']
                elif row['subset'].endswith('Gene=IN'):
                    splcount['IN'] = row['value']

        for cat in CATEGORIES:
            for gene in GENES:
                writer.writerows(prevalence_stat(gene, cat, splcount[gene]))
        # writer.writerows(overall_prevalence_stat())


if __name__ == '__main__':
    main()
