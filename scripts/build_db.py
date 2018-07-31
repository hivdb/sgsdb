#! /usr/bin/env python
"""
Build database files

Please do not call this script directly, use `make build`
"""

import os
import sys
import csv
import json

from collections import OrderedDict


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


def main():
    facttable, fasta, sierra_report, outputdir = sys.argv[1:]
    with open(facttable) as fp1, open(sierra_report) as fp2:
        sequences = list(csv.DictReader(fp1))
        sequence_reports = json.load(fp2)
    sequence_reports = {
        sr['inputSequence']['header'].split('.', 1)[0]: sr
        for sr in sequence_reports
    }
    result_sequences = []
    for seq in sequences:
        seq = OrderedDict(seq)
        accs = seq['Accession']
        seqreport = sequence_reports[accs]
        subtype = seqreport['subtypeText'].split(' (', 1)[0]
        seq.update({
            'PR': 0,
            'RT': 0,
            'IN': 0,
            'Subtype': subtype,
        })
        for gene_seq in seqreport['alignedGeneSequences']:
            gene = gene_seq['gene']['name']
            seq[gene] = 1
        result_sequences.append(seq)
    result_data = {
        'sequences': result_sequences,
        'references': []
    }
    factjson = os.path.join(outputdir, 'meta.json')
    with open(factjson, 'w') as outfp:
        json.dump(result_data, outfp)

    for header, seq in fasta_reader(fasta):
        accs = header.split('.', 1)[0]
        single = os.path.join(outputdir, 'sequences', '{}.json'.format(accs))
        with open(single, 'w') as outfas:
            json.dump({'header': header, 'sequence': seq}, outfas)


if __name__ == '__main__':
    main()
