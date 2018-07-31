sierra:
	@sierrapy fasta -q scripts/query.gql -o local/SGS.sequences.json local/SGS.sequences.fas

build:
	@rm -rf data/upload
	@mkdir -p data/upload/sequences
	@python scripts/build_db.py data/SGS.sequences.fact.csv local/SGS.sequences.fas local/SGS.sequences.json data/upload
