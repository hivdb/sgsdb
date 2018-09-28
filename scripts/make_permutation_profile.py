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
            if row['name'] == '# Patients':
                subset = row['subset']
                if subset == 'Gene=PR':
                    config['PRNumPatients'] = int(row['value'])
                elif subset == 'Gene=RT':
                    config['RTNumPatients'] = int(row['value'])
                elif subset == 'Gene=IN':
                    config['INNumPatients'] = int(row['value'])
                elif subset == 'Subtype=B':
                    config['SubtypeBRatio'] = float(row['percent'][:-1]) / 100
                elif subset == 'Subtype=C':
                    config['SubtypeCRatio'] = float(row['percent'][:-1]) / 100
                elif subset == 'Subtype=Other':
                    config['SubtypeOtherRatio'] = \
                        float(row['percent'][:-1]) / 100
                elif subset == 'Rx=ART':
                    config['RxARTRatio'] = float(row['percent'][:-1]) / 100
                elif subset == 'Rx=None':
                    config['RxNaiveRatio'] = float(row['percent'][:-1]) / 100
    with open(PROFILE_PATH, 'w') as fp:
        json.dump(config, fp, indent=2)


if __name__ == '__main__':
    main()
