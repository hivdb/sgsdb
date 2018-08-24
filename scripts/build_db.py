#! /usr/bin/env python
"""
Build database files

Please do not call this script directly, use `make build`
"""

import os
import sys
import csv
import json
from datetime import datetime

import requests
from collections import OrderedDict

DATE_1900 = datetime.strptime('1900-01-01', '%Y-%m-%d')
ESUMMARY_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'


def fasta_reader(filename):
    with open(filename) as fp:
        header = None
        seq = []
        for line in fp:
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                if seq:
                    yield header, ''.join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            yield header, ''.join(seq)


def retrieve_references(pubmed_ids):
    resp = requests.post(
        ESUMMARY_URL,
        data={
            'db': 'pubmed',
            'id': ','.join(pubmed_ids),
            'retmode': 'json'
        }
    )
    data = resp.json()
    references = OrderedDict()
    for uid in sorted(data['result']['uids']):
        ref = data['result'][uid]
        references[uid] = OrderedDict({
            'MedlineID': uid,
            'PubYear': ref.get('pubdate', '').split(' ')[0],
            'firstAuthor': ref['authors'][0]['name'],
            'Title': ref['title'],
            'Journal': ref['source']
        })
    return references


def list_subtypes(seqs):
    subtypes = set()
    for seq in seqs:
        subtypes.add(seq['Subtype'])
    return sorted(subtypes)


def list_sources(seqs):
    sources = set()
    for seq in seqs:
        sources.add(seq['Source'])
    return sorted(sources)


def main():
    facttable, fasta, sierra_report, outputdir = sys.argv[1:]
    with open(facttable) as fp1, open(sierra_report) as fp2:
        if fp1.read(1) != '\ufeff':
            fp1.seek(0)
        sequences = list(csv.DictReader(fp1))
        sequence_reports = json.load(fp2)
    sequence_reports = {
        sr['inputSequence']['header'].split('.', 1)[0]: sr
        for sr in sequence_reports
    }
    result_sequences = []
    pubmed_ids = set()

    min_days_offsets = {}
    for seq in sequences:
        pt = seq['PtIdentifier']
        offset = \
            datetime.strptime(seq['CollectionDate'], '%Y-%m-%d') - DATE_1900
        offset = offset.days
        min_days_offsets[pt] = \
            min(min_days_offsets.get(pt, 0xffffffff), offset)

    for seq in sequences:
        pt = seq['PtIdentifier']
        offset = \
            datetime.strptime(seq['CollectionDate'], '%Y-%m-%d') - DATE_1900
        offset = offset.days
        weeks = int('{:.0f}'.format((offset - min_days_offsets[pt]) / 7))
        seq = OrderedDict(seq)
        pubmed_ids.add(seq['MedlineID'])
        accs = seq['Accession']
        seqreport = sequence_reports[accs]
        subtype = seqreport['subtypeText'].split(' (', 1)[0]
        seq.update({
            'Weeks': weeks,
            'PR': 0,
            'RT': 0,
            'IN': 0,
            'Subtype': subtype,
        })
        for gene_seq in seqreport['alignedGeneSequences']:
            gene = gene_seq['gene']['name']
            seq[gene] = 1
        result_sequences.append(seq)
    result_data = OrderedDict({
        'sequences': result_sequences,
        'references': retrieve_references(pubmed_ids),
        'subtypes': list_subtypes(result_sequences),
        'sources': list_sources(result_sequences),
    })
    factjson = os.path.join(outputdir, 'meta.json')
    with open(factjson, 'w') as outfp:
        json.dump(result_data, outfp)

    # We don't need this anymore
    # for header, seq in fasta_reader(fasta):
    #     accs = header.split('.', 1)[0]
    #     single = os.path.join(outputdir, 'sequences', '{}.json'.format(accs))
    #     with open(single, 'w') as outfas:
    #         json.dump({'header': header, 'sequence': seq}, outfas)


if __name__ == '__main__':
    main()
