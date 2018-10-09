#! /bin/sh

set -e

pushd `dirname $0` > /dev/null

for gene in PR RT IN; do
    # time python permutation_test.py $gene $1 > ../local/permut.$gene.$1.txt
    time python ordinary_permutation_test.py $gene $1 > ../local/permut.new/permut.$gene.$1.o.txt
done
