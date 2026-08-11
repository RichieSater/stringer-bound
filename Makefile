PYTHON := uv run --frozen python
PY_DIR := supporting-materials/computations/python
TEST_DIR := supporting-materials/computations/tests
CERT_DIR := supporting-materials/computations/certificates
PAPER_DIR := supporting-materials/paper

.PHONY: sync test proof-check all-n-reduction-check dirichlet-poissonization-check poisson-band-certificate-check one-cap-certificate-check n3-formula-check n3-certificate-check n4-structure-check n4-certificate-check n5-structure-check n5-certificate-check certificate-summary-check claim-manifest-check paper reproduce

sync:
	uv sync --frozen

test:
	$(PYTHON) -m unittest discover -s $(TEST_DIR) -v

proof-check:
	$(PYTHON) $(PY_DIR)/n2_proof_check.py

all-n-reduction-check:
	$(PYTHON) $(PY_DIR)/all_n_poisson_reductions.py

dirichlet-poissonization-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/dirichlet_poissonization.py --out $$tmp; \
	cmp $(CERT_DIR)/dirichlet-poissonization-certificate.json $$tmp; \
	rm -f $$tmp

poisson-band-certificate-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/poisson_band_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/poisson-simultaneous-band-certificate.json $$tmp; \
	rm -f $$tmp

one-cap-certificate-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/one_cap_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/one-cap-certificate.json $$tmp; \
	rm -f $$tmp

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

n4-structure-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/derive_n4_bernstein_structure.py --out $$tmp; \
	cmp $(CERT_DIR)/n4-gaffke-bernstein-structure.json $$tmp; \
	rm -f $$tmp

n4-certificate-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/n4_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n4-gaffke-certificate.json $$tmp; \
	rm -f $$tmp

n5-structure-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/derive_n5_bernstein_structure.py --out $$tmp; \
	cmp $(CERT_DIR)/n5-gaffke-bernstein-structure.json $$tmp; \
	rm -f $$tmp

n5-certificate-check:
	@tmp=$$(mktemp); \
	$(PYTHON) $(PY_DIR)/n5_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n5-gaffke-certificate.json $$tmp; \
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

reproduce: sync test proof-check all-n-reduction-check dirichlet-poissonization-check poisson-band-certificate-check one-cap-certificate-check n3-formula-check n3-certificate-check n4-structure-check n4-certificate-check n5-structure-check n5-certificate-check certificate-summary-check claim-manifest-check paper
	@echo "Core proofs, all-n reductions, Dirichlet-Poissonization research certificate, Poisson simultaneous-band and one-cap certificates, n=3 through n=5 Bernstein certificates, counterexample summaries, claim links, tests, and manuscript build passed."
