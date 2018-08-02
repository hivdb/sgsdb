## Steps for updating the SGS database

1. Install `SierraPy` and `NumPy`.
2. Update `data/SGS.sequences.fact.csv` spreadsheet with latest data.
3. Run command `make fasta`; wait until the command finished.
4. Run command `make sierra`; wait until the command finished.
5. Run command `make build stat`;  wait until the command finished.
6. Run command `git add data/upload`, commit and push.
