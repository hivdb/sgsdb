## Steps for updating the SGS database

1. Install `SierraPy`.
2. Update `data/SGS.sequences.fact.csv` spreadsheet with latest data.
3. Create a plain text file contains only the Accession IDs (one per each line).
4. Upload the step 2 file to [Batch Entrz](https://www.ncbi.nlm.nih.gov/sites/batchentrez) and click "retrieve".
5. Select and download the FASTA file.
6. Move the step 4 FASTA file to `local/SGS.sequences.fas`.
7. Run command `make sierra`; wait until the command finished.
8. Run command `make build`;  wait until the command finished.
9. Run command `git add data/upload`, commit and push.
