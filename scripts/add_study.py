#! /usr/bin/env python
from __future__ import print_function

import re
import os
import sys
import csv
import json
from io import BytesIO
from collections import namedtuple

import requests
import sierrapy
from tqdm import tqdm
from Bio import Entrez

Entrez.email = 'hivdbteam@stanford.edu'

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCALDIR = os.path.join(BASEDIR, 'local', 'new_studies')
SIERRAPY_QUERY = os.path.join(BASEDIR, 'scripts', 'query.gql')
LANL_TAG_START = '##HIVDataBaseData-START##'
LANL_TAG_END = '##HIVDataBaseData-END##'

os.makedirs(LOCALDIR, exist_ok=True)


SequenceTuple = namedtuple(
    'SequenceTuple',
    ['accession', 'header', 'sequence',
     'pmid', 'isolate_name', 'isolate_date',
     'source', 'patient', 'extra']
)


def get_accs(pmid):
    resp = requests.get(
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi',
        {'dbfrom': 'pubmed',
         'db': 'nuccore',
         'id': pmid,
         'retmode': 'json'})
    data = json.load(BytesIO(resp.content))
    if 'error' in data:
        return []
    linkset = data['linksets'][0]
    if 'linksetdbs' not in linkset:
        return []
    return linkset['linksetdbs'][0]['links']


def pop_from_qualifier(quals, *names):
    results = []
    for qual in quals:
        if qual['GBQualifier_name'] in names:
            value = qual['GBQualifier_value'].strip()
            results.append(value)
    return results[-1] if results else None


def parse_lanl_data(comment):
    if LANL_TAG_END not in comment:
        return {}
    comment = comment.split(LANL_TAG_END, 1)[0].split(LANL_TAG_START)[1]
    comment = comment.strip(' ;').split(' ; ')
    r = {}
    for n in comment:
        k, v = n.split(' :: ')
        r[k.lower()] = v
    return r


def get_sequences(pmid, accs, step=24):
    print('Fetching sequences from GenBank ...')
    pbar = tqdm(total=len(accs))
    for offset in range(0, len(accs), step):
        partial = accs[offset:offset + step]
        # Entrez.efetch is too slow because it doesn't support gzip
        resp = requests.post(
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
            data={
                'db': 'nuccore',
                'id': ','.join(str(a) for a in partial),
                'retmode': 'xml'
            })
        data = Entrez.read(BytesIO(resp.content))
        for seqdict in data:
            isolate_name = isolate_date = source = patient = None
            extra = {}
            for feature in seqdict['GBSeq_feature-table']:
                if feature['GBFeature_key'] == 'source':
                    quals = feature['GBFeature_quals']
                    isolate_name = pop_from_qualifier(
                        quals, 'isolate', 'strain')
                    isolate_date = pop_from_qualifier(quals, 'collection_date')
                    source = pop_from_qualifier(quals, 'isolation_source')
            if 'GBSeq_comment' in seqdict:
                lanldata = parse_lanl_data(seqdict['GBSeq_comment'])
                source = lanldata.pop('sample tissue', source)
                patient = lanldata.pop('patient code')
                extra = lanldata
            yield SequenceTuple(
                seqdict['GBSeq_primary-accession'],  # accession
                seqdict['GBSeq_accession-version'] +
                ' ' + seqdict['GBSeq_definition'],   # header
                seqdict['GBSeq_sequence'].upper(),   # sequence
                pmid,
                isolate_name,
                isolate_date,
                source,
                patient,
                extra,
            )
        pbar.update(step)


def get_sierra_result(sequences):
    seqs = [{'header': s.header, 'sequence': s.sequence}
            for s in sequences]
    client = sierrapy.SierraClient('https://hivdb.stanford.edu/graphql')
    client.toggle_progress()
    with open(SIERRAPY_QUERY) as fp:
        query = fp.read()
    return client.sequence_analysis(seqs, query)


def write_sequences_fact(sequences, sierra_result, fp):
    accs = set()
    for one in sierra_result:
        if not one['alignedGeneSequences']:
            continue
        accs.add(one['inputSequence']['header'].split('.', 1)[0])

    extrakeys = set()
    for seq in sequences:
        if seq.accession not in accs:
            continue
        extrakeys |= seq.extra.keys()
    extrakeys = sorted(extrakeys)

    fields = [
        'MedlineID', 'Accession', 'CollectionDate',
        'IsolateName', 'Patient', 'Source', 'Definition',
        *['extra.{}'.format(k) for k in extrakeys]
    ]
    writer = csv.writer(fp)
    writer.writerow(fields)
    for seq in sequences:
        if seq.accession not in accs:
            continue
        writer.writerow([
            seq.pmid, seq.accession, seq.isolate_date or '',
            seq.isolate_name or '', seq.patient or '',
            seq.source or '', seq.header,
            *[seq.extra[k] for k in extrakeys]
        ])
    if accs:
        print('{} non-pol sequences are removed from the final result.'
              .format(len(sequences) - len(accs)), file=sys.stderr)


def write_sierra_result(sierra_result, fp):
    data = [o for o in sierra_result if o['alignedGeneSequences']]
    json.dump(data, fp, indent=2)


def single(pmid, accs=None):
    if not re.match(r'^\d+$', pmid):
        print('PMID must be a number (received {!r}'
              .format(pmid), file=sys.stderr)
        exit(2)

    # fetch accessions
    if not accs:
        print('Fetching GenBank accession number from Entrez ...',
              file=sys.stderr, end='\r')
        accs = get_accs(pmid)
        if accs:
            print('Fetching GenBank accession number from Entrez:',
                  'found {}.'.format(len(accs)), file=sys.stderr)
        else:
            print('Fetching GenBank accession number from Entrez:',
                  'not found.', file=sys.stderr)
            exit(3)

    # fetch sequences
    sequences = list(get_sequences(pmid, accs))

    # run sierra
    sierra_result = get_sierra_result(sequences)

    # save sequences fact
    fact_table = os.path.join(LOCALDIR, '{}.sequences.fact.csv'.format(pmid))
    with open(fact_table, 'w') as fp:
        write_sequences_fact(sequences, sierra_result, fp)

    # save sierra json
    sierra_file = os.path.join(LOCALDIR, '{}.sequences.json'.format(pmid))
    with open(sierra_file, 'w') as fp:
        write_sierra_result(sierra_result, fp)
    print('- {}'.format(sierra_file), file=sys.stderr)
    print('- {}'.format(fact_table), file=sys.stderr)


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ('multiple', 'single'):
        print('Usage: {} <multiple|single> ...'
              .format(sys.argv[0]), file=sys.stderr)
        exit(127)
    method = sys.argv[1]
    if method == 'multiple':
        if len(sys.argv) < 3:
            print('Usage: {} multiple <PMID1> [PMID2, PMID3, ...]'
                  .format(sys.argv[0]), file=sys.stderr)
            exit(126)
        pmids = sys.argv[2:]
        for pmid in pmids:
            single(pmid)
    else:  # single
        if len(sys.argv) < 3:
            print('Usage: {} single <PMID> [ACCESSION1, ACCESSION2, ...]'
                  .format(sys.argv[0]), file=sys.stderr)
            exit(125)
        (pmid, *accs) = sys.argv[2:]
        single(pmid, accs)


if __name__ == '__main__':
    main()
