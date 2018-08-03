#! /usr/bin/env python

import os
import json
import csv

from itertools import groupby
from numpy import percentile
from collections import defaultdict, Counter

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)

FACTSHEET = os.path.join(BASEDIR, 'data', 'SGS.sequences.fact.csv')
SIEERAREPORT = os.path.join(BASEDIR, 'local', 'SGS.sequences.json')
OUTPUT = os.path.join(BASEDIR, 'data', 'SGS.study-summary.csv')

OUTPUT_HEADER = [
    'Study', '# Patients', '# Sequences', 'Subtypes',
    '# Patients per Subtype and Gene',
    '# Sequences per Subtype and Gene', 'Mixtures',
    'PR GeneSize', 'RT GeneSize', 'IN GeneSize',
    'PR FirstAA', 'PR LastAA',
    'RT FirstAA', 'RT LastAA',
    'IN FirstAA', 'IN LastAA',
]


def sorted_groupby(iterable, keyfunc=None):
    iterable = sorted(iterable, key=keyfunc)
    return groupby(iterable, keyfunc)


def rangetext(arr, float_format='.0f'):
    if not arr:
        return 'N/A'
    left, med, right = percentile(arr, (10, 50, 90))
    tpl = '{{:{0}}} ({{:{0}}} - {{:{0}}})'.format(float_format)
    return tpl.format(med, left, right)


def counter2text(counter, order_by_name=False):
    result = []
    if order_by_name:
        items = sorted(counter.items())
    else:
        items = counter.most_common()
    for name, num in items:
        result.append('{} ({})'.format(name, num))
    return ', '.join(result)


def sequence_set_summary(sequences, sequence_reports):
    genesizes = defaultdict(list)
    first_aas = defaultdict(list)
    last_aas = defaultdict(list)
    sequences = list(sequences)
    subtypes = Counter()
    all_mixtures = Counter()
    pt_per_subtype_gene = defaultdict(set)
    seq_per_subtype_gene = Counter()
    sequences = [s for s in sequences if s['Source'] == 'Plasma']
    for seq in sequences:
        n_mixtures = 0
        report = sequence_reports[seq['Accession']]
        subtype = report['subtypeText'].split(' (', 1)[0]
        subtypes[subtype] += 1
        for geneseq in report['alignedGeneSequences']:
            gene = geneseq['gene']['name']

            first_aa = int(geneseq['firstAA'])
            last_aa = int(geneseq['lastAA'])
            genesize = last_aa - first_aa + 1
            genesizes[gene].append(genesize)
            first_aas[gene].append(first_aa)
            last_aas[gene].append(last_aa)
            for mut in geneseq['mutations']:
                aas = mut['AAs']
                if '-' in aas:
                    continue
                aas = aas.split('_', 1)[0]
                n_mixtures += len(aas) > 1
            (pt_per_subtype_gene['{} {}'.format(subtype, gene)]
             .add(seq['PtIdentifier']))
            (pt_per_subtype_gene['All {}'.format(gene)]
             .add(seq['PtIdentifier']))
            seq_per_subtype_gene['{} {}'.format(subtype, gene)] += 1
            seq_per_subtype_gene['All {}'.format(gene)] += 1
        all_mixtures[n_mixtures] += 1
    pt_per_subtype_gene = Counter(
        {k: len(v) for k, v in pt_per_subtype_gene.items()})
    return {
        '# Patients': len(set(s['PtIdentifier'] for s in sequences)),
        '# Sequences': len(set(s['Accession'] for s in sequences)),
        'Subtypes': counter2text(subtypes),
        '# Patients per Subtype and Gene': counter2text(pt_per_subtype_gene),
        '# Sequences per Subtype and Gene': counter2text(seq_per_subtype_gene),
        'Mixtures': counter2text(all_mixtures, order_by_name=True),
        'PR GeneSize': rangetext(genesizes['PR']),
        'PR FirstAA': rangetext(first_aas['PR']),
        'PR LastAA': rangetext(last_aas['PR']),
        'RT GeneSize': rangetext(genesizes['RT']),
        'RT FirstAA': rangetext(first_aas['RT']),
        'RT LastAA': rangetext(last_aas['RT']),
        'IN GeneSize': rangetext(genesizes['IN']),
        'IN FirstAA': rangetext(first_aas['IN']),
        'IN LastAA': rangetext(last_aas['IN']),
    }


def main():
    with open(FACTSHEET) as fp:
        if fp.read(1) != '\ufeff':
            fp.seek(0)
        sequences = list(csv.DictReader(fp))
    with open(SIEERAREPORT) as fp:
        sequence_reports = json.load(fp)
    sequence_reports = {s['inputSequence']['header'].split('.', 1)[0]: s
                        for s in sequence_reports}
    with open(OUTPUT, 'w') as fp:
        writer = csv.DictWriter(fp, OUTPUT_HEADER)
        writer.writeheader()
        overalls = sequence_set_summary(sequences, sequence_reports)
        overalls['Study'] = 'Overall'
        for medid, study_sequences in \
                sorted_groupby(sequences, lambda s: s['MedlineID']):
            summary = sequence_set_summary(study_sequences, sequence_reports)
            summary['Study'] = medid
            writer.writerow(summary)
        writer.writerow(overalls)


if __name__ == '__main__':
    main()
