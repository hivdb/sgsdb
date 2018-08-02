#! /usr/bin/env python

import os
import csv
import json

from decimal import Decimal
from collections import OrderedDict, Counter, defaultdict

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)

FACTSHEET = os.path.join(BASEDIR, 'SGS.sequences.fact.csv')
SIEERAREPORT = os.path.join(BASEDIR, 'SGS.sequences.json')
OUTPUTS = {
    'PR': os.path.join(BASEDIR, 'prevalenceData', 'SGS.PRprevalence.csv'),
    'RT': os.path.join(BASEDIR, 'prevalenceData', 'SGS.RTprevalence.csv'),
    'IN': os.path.join(BASEDIR, 'prevalenceData', 'SGS.INprevalence.csv'),
}

PREC3 = Decimal('1.000')

CONSENSUS = {
    'PR': (
        'PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGH'
        'KAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'
    ),
    'RT': (
        'PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDST'
        'KWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNE'
        'TPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQ'
        'HLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYAGIKV'
        'KQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEP'
        'FKNLKTGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEAWWTEYWQATWIPEW'
        'EFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTDRGRQKVVSLTDTTNQKTELQAIH'
        'LALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLV'
        'SAGIRKVL'
    ),
    'IN': (
        'FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLE'
        'GKIILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQE'
        'FGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATD'
        'IQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMA'
        'GDDCVASRQDED'
    ),
}


def aggregate_aa_prevalence(gene, sequences, seqsfact, subtype=None):
    result = OrderedDict()
    total = Counter()
    totalpt = defaultdict(set)
    resultpt = defaultdict(set)
    all_aas = list('ACDEFGHIKLMNPQRSTVWY') + ['del', 'X', 'ins', '*']
    consensus = CONSENSUS[gene]
    genesize = len(consensus)

    for pos in range(1, genesize + 1):
        for aa in all_aas:
            result[(pos, aa)] = 0

    for seq in sequences:
        accn = seq['inputSequence']['header'].split('.', 1)[0]
        seqfact = seqsfact[accn]
        seq_subtype = seq['subtypeText'].split(' (', 1)[0]
        if subtype is not None and subtype != seq_subtype:
            continue
        for gseq in seq['alignedGeneSequences']:
            if gseq['gene']['name'] == gene:
                break
        else:
            continue
        muts = {m['position']: m for m in gseq['mutations']}
        for pos in range(1, genesize + 1):
            cons = consensus[pos - 1]
            if pos in muts:
                mut = muts[pos]
                if mut['isInsertion']:
                    aa = 'ins'
                elif mut['isDeletion']:
                    aa = 'del'
                else:
                    aa = mut['AAs'].replace('cons', '')
                    if len(aa) > 1:
                        aa = 'X'
            else:
                aa = cons
            result[(pos, aa)] += 1
            total[pos] += 1
            resultpt[(pos, aa)].add(seqfact['PtIdentifier'])
            totalpt[pos].add(seqfact['PtIdentifier'])
    return OrderedDict(((pos, aa), {
        'Gene': gene,
        'Subtype': subtype,
        'Pos': pos,
        'AA': aa,
        'Pcnt': Decimal(count /
                        (total[pos] or 0.001) * 100).quantize(PREC3),
        'Count': count,
        'PosTotal': total[pos],
        'PatientCount': len(resultpt[(pos, aa)]),
        'PatientPosTotal': len(totalpt[pos]),
    }) for (pos, aa), count in result.items())


def main():
    major_subtypes = [None, 'B', 'C', 'D', 'A']
    header = ['Gene', 'Subtype', 'Pos', 'AA', 'Pcnt', 'Count',
              'PosTotal', 'PatientCount', 'PatientPosTotal']
    with open(SIEERAREPORT) as fp:
        sequences = json.load(fp)
    with open(FACTSHEET) as fp:
        if fp.read(1) != '\ufeff':
            fp.seek(0)
        seqsfact = csv.DictReader(fp)
        seqsfact = {s['Accession']: s for s in seqsfact}
    for gene in ('PR', 'RT', 'IN'):
        all_prevalence = []
        for subtype in major_subtypes:
            prevs = aggregate_aa_prevalence(
                gene, sequences, seqsfact, subtype).values()
            all_prevalence.extend(prevs)
        with open(OUTPUTS[gene], 'w') as fp:
            writer = csv.DictWriter(fp, header)
            writer.writeheader()
            writer.writerows(all_prevalence)


if __name__ == '__main__':
    main()
