PYTHON := uv run --frozen python
PY_DIR := supporting-materials/computations/python
TEST_DIR := supporting-materials/computations/tests
CERT_DIR := supporting-materials/computations/certificates
PAPER_DIR := supporting-materials/paper

.PHONY: sync test proof-check n3-formula-check n3-certificate-check certificate-summary-check claim-manifest-check paper reproduce

sync:
	uv sync --frozen

test:
	$(PYTHON) -m unittest discover -s $(TEST_DIR) -v

proof-check:
	$(PYTHON) $(PY_DIR)/n2_proof_check.py

n3-formula-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/derive_n3_bernstein_formulas.py --out $$tmp; \
	cmp $(CERT_DIR)/n3-gaffke-bernstein-formulas.json $$tmp; \
	rm -f $$tmp

n3-certificate-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/n3_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n3-gaffke-certificate.json $$tmp; \
	rm -f $$tmp

certificate-summary-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/summarize_certificates.py --out $$tmp >/dev/null; \
	cmp $(CERT_DIR)/certificate-summary.json $$tmp; \
	rm -f $$tmp

claim-manifest-check:
	$(PYTHON) $(PY_DIR)/validate_claim_manifest.py

paper:
	cd $(PAPER_DIR) && tectonic -X compile stringer.tex --keep-logs --keep-intermediates
	@! grep -Eq "Warning|Overfull|Underfull|undefined|multiply defined" $(PAPER_DIR)/stringer.log

reproduce: sync test proof-check n3-formula-check n3-certificate-check certificate-summary-check claim-manifest-check paper
	@echo "Core proofs, n=3 Bernstein certificate, counterexample summaries, claim links, tests, and manuscript build passed."
