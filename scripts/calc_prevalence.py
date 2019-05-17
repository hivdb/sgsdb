#! /usr/bin/env python

import os
import csv

from decimal import Decimal
from collections import OrderedDict, Counter, defaultdict

from common import load_sequences, apobec_mutation_map, BASEDIR, PREC3

OUTPUTS = {
    'PR': os.path.join(BASEDIR, 'data', 'prevalence', 'SGS.PRprevalence.csv'),
    'RT': os.path.join(BASEDIR, 'data', 'prevalence', 'SGS.RTprevalence.csv'),
    'IN': os.path.join(BASEDIR, 'data', 'prevalence', 'SGS.INprevalence.csv'),
}

APM = apobec_mutation_map()

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

NA_CONSENSUS = {
    'PR': (
        'CCTCAGATCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTA'
        'TTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATG'
        'ATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACAT'
        'AAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAG'
        'ATTGGTTGCACTTTAAATTTT'
    ),
    'RT': (
        'CCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAA'
        'TGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGGAAA'
        'ATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACT'
        'AAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTA'
        'GGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATAT'
        'TTTTCAGTTCCCTTAGATAAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAG'
        'ACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAA'
        'AGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATG'
        'GATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAA'
        'CATCTGTTGAGGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATG'
        'GGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACT'
        'GTCAATGACATACAGAAGTTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTATGCAGGGATTAAAGTA'
        'AAGCAATTATGTAAACTCCTTAGGGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCA'
        'GAGCTAGAACTGGCAGAAAACAGGGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCA'
        'AAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCA'
        'TTTAAAAATCTGAAAACAGGAAAGTATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTA'
        'ACAGAGGCAGTGCAAAAAATAGCCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTA'
        'CCCATACAAAAAGAAACATGGGAAGCATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGG'
        'GAGTTTGTCAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCA'
        'GAAACTTTCTATGTAGATGGGGCAGCTAATAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTGAC'
        'AGAGGAAGACAAAAAGTTGTCTCCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTCAT'
        'CTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATT'
        'CAAGCACAACCAGATAAAAGTGAATCAGAGTTAGTCAGTCAAATAATAGAGCAGTTAATAAAAAAGGAA'
        'AAGGTCTACCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTC'
        'AGTGCTGGAATCAGGAAAGTACTA'
    ),
    'IN': (
        'TTTTTAGATGGAATAGATAAGGCCCAAGAAGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCT'
        'AGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAA'
        'GGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAA'
        'GGAAAAATTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAG'
        'ACAGGGCAGGAAACAGCATACTTTCTCTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACA'
        'GACAATGGCAGCAATTTCACCAGTACTACGGTTAAGGCCGCCTGTTGGTGGGCAGGGATCAAGCAGGAA'
        'TTTGGCATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATT'
        'ATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAAT'
        'TTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGAC'
        'ATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGC'
        'AGAGATCCACTTTGGAAAGGACCAGCAAAGCTTCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGAT'
        'AATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCA'
        'GGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT'
    )
}


def displaycodons(codoncounter):
    tpl = '{} ({})'
    if codoncounter and codoncounter.most_common(1)[0][1] == 1:
        tpl = '{}'
    return ', '.join(tpl.format(*c) for c in codoncounter.most_common())


def othercodons(ptids, codons, exclude_codons):
    result = Counter()
    for ptid, ptcodons in codons.items():
        if ptid in ptids:
            result += ptcodons
    for codon in exclude_codons:
        result.pop(codon)
    return displaycodons(result)


def aggregate_aa_prevalence(gene, sequences,
                            category, category_func):
    result = OrderedDict()
    total = Counter()
    totalpt = defaultdict(set)
    resultpt = defaultdict(set)
    totalspl = defaultdict(set)
    resultspl = defaultdict(set)
    all_aas = list('ACDEFGHIKLMNPQRSTVWY-_X*')
    consensus = CONSENSUS[gene]
    genesize = len(consensus)
    codons = defaultdict(Counter)
    shortcodons = defaultdict(Counter)
    codonspt = defaultdict(lambda: defaultdict(Counter))
    shortcodonspt = defaultdict(lambda: defaultdict(Counter))

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
        aligned_nas = gseq['alignedNAs']
        first_aa = gseq['firstAA']
        last_aa = gseq['lastAA']
        for pos in range(1, genesize + 1):
            if pos < first_aa or pos > last_aa:
                continue
            cons = consensus[pos - 1]
            if pos in muts:
                mut = muts[pos]
                if mut['isInsertion']:
                    aa = '_'
                elif mut['isDeletion']:
                    aa = '-'
                else:
                    aa = mut['AAs']
                    if len(aa) > 1:
                        aa = 'X'
            else:
                aa = cons
            result[(pos, aa)] += 1
            relpos = pos - first_aa + 1
            total[pos] += 1
            codon = aligned_nas[max(0, (relpos - 2) * 3):(relpos + 1) * 3]
            shortcodon = aligned_nas[relpos * 3 - 3:relpos * 3]
            ptid = seqfact['PtIdentifier']
            tp = seqfact['CollectionDate']
            codons[(pos, aa)][codon] += 1
            shortcodons[(pos, aa)][shortcodon] += 1
            codonspt[pos][ptid][codon] += 1
            shortcodonspt[pos][ptid][shortcodon] += 1
            resultpt[(pos, aa)].add(ptid)
            totalpt[pos].add(ptid)
            resultspl[(pos, aa)].add((ptid, tp))
            totalspl[pos].add((ptid, tp))
    return OrderedDict(((pos, aa), {
        'Gene': gene,
        'Category': category,
        'Pos': pos,
        'AA': aa,
        # 'NACons': NA_CONSENSUS[gene][max(0, (pos - 2) * 3):(pos + 1) * 3],
        'FromCodons': othercodons(
            resultpt[(pos, aa)], shortcodonspt[pos], shortcodons[(pos, aa)]
        ),
        'ToCodons': displaycodons(shortcodons[(pos, aa)]),
        'FromCodonsCtx': othercodons(
            resultpt[(pos, aa)], codonspt[pos], codons[(pos, aa)]
        ),
        'ToCodonsCtx': displaycodons(codons[(pos, aa)]),
        'Pcnt': Decimal(count /
                        (total[pos] or 0.001) * 100).quantize(PREC3),
        'Count': count,
        'PosTotal': total[pos],
        'PatientCount': len(resultpt[(pos, aa)]),
        'PatientPosTotal': len(totalpt[pos]),
        'SampleCount': len(resultspl[(pos, aa)]),
        'SamplePosTotal': len(totalspl[pos]),
        'IsAPOBEC': (gene, pos, aa) in APM,
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
    header = ['Gene', 'Category', 'Pos', 'AA', 'FromCodons', 'ToCodons',
              'FromCodonsCtx', 'ToCodonsCtx', 'Pcnt', 'Count', 'PosTotal',
              'PatientCount', 'PatientPosTotal', 'SampleCount',
              'SamplePosTotal', 'IsAPOBEC']
    sequences = load_sequences(filtered=True)
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
