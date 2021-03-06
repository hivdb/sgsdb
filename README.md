## Steps for updating the SGS database

1. Install dependencies with `pipenv install`. Python 3.9 is required.
2. Update `data/SGS.sequences.fact.csv` spreadsheet with latest data. Follow [Steps
   for adding studies](#steps-for-adding-studies) for creating intermediate files.
3. Run command `make fasta`; wait until the command finished.
4. Run command `make sierra`; wait until the command finished.
5. Run command `make build stat`;  wait until the command finished.
6. Run command `git add data/upload`, commit and push.

## Steps for adding studies

```
pipenv run python scripts/add_study.py multiple <PMID1> <PMID2> <PMID3> ...

# or

pipenv run python scripts/add_study.py single <PMID> [ACCESSION1, ACCESSION2, ...]

```

Alternative way:

1. Find Genbank IDs for this study.
2. In Mangabey data entry program, check "No reference" checkbox in reference
   entry page then "Continue".
3. Click "Nucleotide Sequences".
4. Click "Add Genbank Sequence(s)" then type the Genbank IDs. Wait until all
   sequences loaded.
5. "Download" the output file.
6. Merge the downloaded TSV file with `data/SGS.sequences.fact.csv`. Remember
   to delete all non-pol sequences.

### Fields

- MedlineID
- Accession
- CollectionDate: Format YYYY-MM-DD
- Source: Plasma, PBMC, etc
- PtIdentifier: Format '{InternalRefID}-{SourcePtID}'
- CometSubtype: Use COMET to determine subtype; can be empty
- Rx: ART or None
- DateAdded: Format YYYY-MM-DD
- _Include: Always False (only used in research)
- _Reservoir: Is this sample collected from virus reservoir?
