#! /usr/bin/env python
import os
import csv
import json

from common import BASEDIR

REPORT_PATH = os.path.join(BASEDIR, 'data', 'report.csv')
PROFILE_PATH = os.path.join(BASEDIR, 'local', 'permutation_profile.json')


def main():
    config = {}
    with open(REPORT_PATH) as fp:
        data = csv.DictReader(fp)
        for row in data:
            if row['name'] == '# Samples (Patient Time Points)':
                subset = row['subset']
                if subset == 'SequencesPerSample>9, Gene=PR':
                    config['PRNumSamples'] = int(row['value'])
                elif subset == 'SequencesPerSample>9, Gene=RT':
                    config['RTNumSamples'] = int(row['value'])
                elif subset == 'SequencesPerSample>9, Gene=IN':
                    config['INNumSamples'] = int(row['value'])
                elif subset == 'SequencesPerSample>9, Subtype=B':
                    config['SubtypeBRatio'] = float(row['percent'][:-1]) / 100
                elif subset == 'SequencesPerSample>9, Subtype=C':
                    config['SubtypeCRatio'] = float(row['percent'][:-1]) / 100
                elif subset == 'SequencesPerSample>9, Subtype=Other':
                    config['SubtypeOtherRatio'] = \
                        float(row['percent'][:-1]) / 100
                elif subset == 'SequencesPerSample>9, Rx=ART':
                    config['RxARTRatio'] = float(row['percent'][:-1]) / 100
                elif subset == 'SequencesPerSample>9, Rx=None':
                    config['RxNaiveRatio'] = float(row['percent'][:-1]) / 100
    with open(PROFILE_PATH, 'w') as fp:
        json.dump(config, fp, indent=2)


if __name__ == '__main__':
    main()
