.PHONY: help test smoke golden

help:
	@echo "make test    - fast pure-logic unit tests (no network, no databases)"
	@echo "make smoke   - end-to-end smoke test (needs the hmm-discovery conda env)"
	@echo "make golden DISCOVERY=<run_dir> FAMILY=<name>  - check a run against tests/golden/<name>.json"

test:
	python3 scripts/run_tests.py

smoke:
	bash run.sh --fasta examples/example_seeds.fasta --smoke --no-controls

golden:
	python3 scripts/golden_test.py --discovery "$(DISCOVERY)" --baseline tests/golden/$(FAMILY).json
