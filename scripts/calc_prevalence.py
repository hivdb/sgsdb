#! /usr/bin/env python

import os
import csv

from decimal import Decimal
from collections import OrderedDict, Counter, defaultdict

from common import load_sequences, BASEDIR, PREC3

OUTPUTS = {
    'PR': os.path.join(BASEDIR, 'data', 'prevalence', 'SGS.PRprevalence.csv'),
    'RT': os.path.join(BASEDIR, 'data', 'prevalence', 'SGS.RTprevalence.csv'),
    'IN': os.path.join(BASEDIR, 'data', 'prevalence', 'SGS.INprevalence.csv'),
}


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


def aggregate_aa_prevalence(gene, sequences,
                            category, category_func):
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

    for seqfact in sequences:
        seq = seqfact['_Sierra']
        seq_subtype = seq['subtypeText'].split(' (', 1)[0]
        if not category_func(seq_subtype, seqfact['Rx']):
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
        'Category': category,
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
    categories = {
        'All': lambda s, rx: True,
        'SubtypeB': lambda s, rx: s == 'B',
        'SubtypeC': lambda s, rx: s == 'C',
        'Non-SubtypeBC': lambda s, rx: s not in ('B', 'C'),
        'ART': lambda s, rx: rx == 'ART',
        'Naive': lambda s, rx: rx == 'None'
    }
    header = ['Gene', 'Category', 'Pos', 'AA', 'Pcnt', 'Count',
              'PosTotal', 'PatientCount', 'PatientPosTotal']
    sequences = load_sequences()
    for gene in ('PR', 'RT', 'IN'):
        all_prevalence = []
        for cat, func in categories.items():
            prevs = aggregate_aa_prevalence(
                gene, sequences, cat, func).values()
            all_prevalence.extend(prevs)
        with open(OUTPUTS[gene], 'w') as fp:
            writer = csv.DictWriter(fp, header)
            writer.writeheader()
            writer.writerows(all_prevalence)


if __name__ == '__main__':
    main()
