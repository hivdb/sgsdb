fasta:
	@rm local/SGS.sequences.fas 2>/dev/null || true
	@scripts/download_fasta.sh

sierra:
	@sierrapy fasta -q scripts/query.gql -o local/SGS.sequences.json local/SGS.sequences.fas
	@zip -j -FS data/SGS.sequences.json.zip local/SGS.sequences.json

build:
	@rm -rf data/upload
	@mkdir -p data/upload
	@python scripts/build_db.py data/SGS.sequences.fact.csv local/SGS.sequences.fas local/SGS.sequences.json data/upload

stat:
	@mkdir -p data/prevalence
	@python scripts/calc_prevalence.py
	@Rscript scripts/comparePrevalence.r
	@python scripts/make_report.py

permutation:
	@python scripts/make_permutation_profile.py
	@scripts/permutation_test.sh 1000

permutation_fig:
	@Rscript scripts/permutationTestGraphOrdinaryMuts.R
