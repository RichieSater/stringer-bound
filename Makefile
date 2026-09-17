PYTHON := uv run --frozen python
PY_DIR := supporting-materials/computations/python
TEST_DIR := supporting-materials/computations/tests
CERT_DIR := supporting-materials/computations/certificates
PAPER_DIR := supporting-materials/paper

.PHONY: sync test public-corpus-check proof-check all-n-reduction-check tetrahedral-vertex-barrier-check five-coordinate-mixed-section-check five-coordinate-block-faces-check six-coordinate-repeated-top-check six-coordinate-unique-top-check six-coordinate-unique-top-l2-l3-check n4-weight-monotonicity-check n5-weight-monotonicity-check terminal-edge-bifurcation-check three-exponential-interior-check srswor-conditioning-check systematic-pps-check systematic-minimax-check systematic-disjoint-dp-check dirichlet-poissonization-check poisson-band-certificate-check poisson-band-calibration-check audit-benchmark-check simplex-cap-interface-check one-cap-all-n-check one-cap-certificate-check two-cap-certificate-check n3-formula-check n3-certificate-check n4-structure-check n4-certificate-check n5-structure-check n5-certificate-check n6-structure-check n6-structure-data-check n6-face-standard-check n6-face-standard-0-check n6-face-standard-1-check n6-face-standard-2-check n6-face-standard-3-check n6-face-standard-4-check n6-face-standard-5-check n6-face-c6-check n6-face-c6-0-check n6-face-c6-1-check n6-face-c6-2-check n6-face-d6-check n6-face-d6-0-check n6-face-d6-1-check n6-face-d6-2-check n6-certificate-check n6-certificate-001-check n6-certificate-005-check n6-certificate-010-check n7-structure-data-check n7-factor-order-check n7-face-a-check n7-face-b-basic-check n7-face-b5-check n7-face-c-basic-check n7-face-c6-check n7-face-c8-check n7-face-d-basic-check n7-face-d6-check n7-face-d9-check n7-certificate-check n7-certificate-001-check n7-certificate-005-check n7-certificate-010-check certificate-summary-check claim-manifest-check paper reproduce

sync:
	uv sync --frozen

test:
	$(PYTHON) -m unittest discover -s $(TEST_DIR) -v

public-corpus-check:
	$(PYTHON) $(PY_DIR)/check_public_corpus.py

proof-check:
	$(PYTHON) $(PY_DIR)/n2_proof_check.py

all-n-reduction-check:
	$(PYTHON) $(PY_DIR)/all_n_poisson_reductions.py

tetrahedral-vertex-barrier-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/tetrahedral_vertex_barrier.py --out $$tmp; \
	cmp $(CERT_DIR)/tetrahedral-vertex-barrier-certificate.json $$tmp

five-coordinate-mixed-section-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/five_coordinate_mixed_section.py --out $$tmp; \
	cmp $(CERT_DIR)/five-coordinate-mixed-section-certificate.json $$tmp

five-coordinate-block-faces-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/five_coordinate_block_faces.py --out $$tmp; \
	cmp $(CERT_DIR)/five-coordinate-block-faces-certificate.json $$tmp

six-coordinate-repeated-top-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/six_coordinate_repeated_top.py --out $$tmp; \
	cmp $(CERT_DIR)/six-coordinate-repeated-top-certificate.json $$tmp

six-coordinate-unique-top-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/six_coordinate_unique_top.py --out $$tmp; \
	cmp $(CERT_DIR)/six-coordinate-unique-top-certificate.json $$tmp

six-coordinate-unique-top-l2-l3-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/six_coordinate_unique_top_l2_l3.py --out $$tmp; \
	cmp $(CERT_DIR)/six-coordinate-unique-top-l2-l3-certificate.json $$tmp

n4-weight-monotonicity-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/n4_weight_monotonicity.py --out $$tmp; \
	cmp $(CERT_DIR)/n4-weight-monotonicity-certificate.json $$tmp

n5-weight-monotonicity-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/n5_weight_monotonicity.py --out $$tmp; \
	cmp $(CERT_DIR)/n5-weight-monotonicity-certificate.json $$tmp

terminal-edge-bifurcation-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/terminal_edge_bifurcation.py --out $$tmp; \
	cmp $(CERT_DIR)/terminal-edge-bifurcation-certificate.json $$tmp

three-exponential-interior-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/three_exponential_interior_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/three-exponential-interior-certificate.json $$tmp

srswor-conditioning-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/srswor_conditioning.py --out $$tmp; \
	cmp $(CERT_DIR)/srswor-conditioning-certificate.json $$tmp

systematic-pps-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/systematic_pps.py --out $$tmp; \
	cmp $(CERT_DIR)/systematic-pps-certificate.json $$tmp

systematic-minimax-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/systematic_minimax.py --out $$tmp; \
	cmp $(CERT_DIR)/systematic-minimax-certificate.json $$tmp

systematic-disjoint-dp-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/systematic_disjoint_dp.py --out $$tmp; \
	cmp $(CERT_DIR)/systematic-disjoint-dp-certificate.json $$tmp

dirichlet-poissonization-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/dirichlet_poissonization.py --out $$tmp; \
	cmp $(CERT_DIR)/dirichlet-poissonization-certificate.json $$tmp

poisson-band-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/poisson_band_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/poisson-simultaneous-band-certificate.json $$tmp

poisson-band-calibration-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/poisson_band_calibration.py --out $$tmp; \
	cmp $(CERT_DIR)/poisson-band-calibration-certificate.json $$tmp

audit-benchmark-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/audit_decision_benchmark.py --out $$tmp; \
	cmp $(CERT_DIR)/audit-decision-benchmark.json $$tmp

simplex-cap-interface-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/simplex_cap_interface.py --out $$tmp; \
	cmp $(CERT_DIR)/simplex-cap-solution-index.json $$tmp

one-cap-all-n-check:
	$(PYTHON) $(PY_DIR)/one_cap_all_n_check.py

one-cap-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/one_cap_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/one-cap-certificate.json $$tmp

two-cap-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/two_cap_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/two-cap-certificate.json $$tmp

n3-formula-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/derive_n3_bernstein_formulas.py --out $$tmp; \
	cmp $(CERT_DIR)/n3-gaffke-bernstein-formulas.json $$tmp

n3-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/n3_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n3-gaffke-certificate.json $$tmp

n4-structure-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/derive_n4_bernstein_structure.py --out $$tmp; \
	cmp $(CERT_DIR)/n4-gaffke-bernstein-structure.json $$tmp

n4-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/n4_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n4-gaffke-certificate.json $$tmp

n5-structure-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/derive_n5_bernstein_structure.py --out $$tmp; \
	cmp $(CERT_DIR)/n5-gaffke-bernstein-structure.json $$tmp

n5-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/n5_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n5-gaffke-certificate.json $$tmp

n6-structure-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py --out $$tmp; \
	cmp $(CERT_DIR)/n6-gaffke-bernstein-structure.json $$tmp

n6-structure-data-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--check-data-against $(CERT_DIR)/n6-gaffke-bernstein-structure.json

n6-face-standard-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group standard

n6-face-standard-0-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group standard-0

n6-face-standard-1-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group standard-1

n6-face-standard-2-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group standard-2

n6-face-standard-3-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group standard-3

n6-face-standard-4-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group standard-4

n6-face-standard-5-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group standard-5

n6-face-c6-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group c6

n6-face-c6-0-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group c6-0

n6-face-c6-1-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group c6-1

n6-face-c6-2-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group c6-2

n6-face-d6-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group d6

n6-face-d6-0-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group d6-0

n6-face-d6-1-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group d6-1

n6-face-d6-2-check:
	$(PYTHON) $(PY_DIR)/derive_n6_bernstein_structure.py \
		--face-proof-group d6-2

n6-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/n6_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n6-gaffke-certificate.json $$tmp

n6-certificate-001-check:
	$(PYTHON) $(PY_DIR)/n6_gaffke_certificate.py --check-alpha 0.01

n6-certificate-005-check:
	$(PYTHON) $(PY_DIR)/n6_gaffke_certificate.py --check-alpha 0.05

n6-certificate-010-check:
	$(PYTHON) $(PY_DIR)/n6_gaffke_certificate.py --check-alpha 0.10

n7-structure-data-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--check-data-against $(CERT_DIR)/n7-gaffke-bernstein-structure.json.gz

n7-factor-order-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--factor-order-check

n7-face-a-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group a

n7-face-b-basic-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group b-basic

n7-face-b5-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group b5

n7-face-c-basic-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group c-basic

n7-face-c6-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group c6

n7-face-c8-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group c8

n7-face-d-basic-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group d-basic

n7-face-d6-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group d6

n7-face-d9-check:
	$(PYTHON) $(PY_DIR)/derive_n7_bernstein_structure.py \
		--face-proof-group d9

n7-certificate-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/n7_gaffke_certificate.py --out $$tmp; \
	cmp $(CERT_DIR)/n7-gaffke-certificate.json $$tmp

n7-certificate-001-check:
	$(PYTHON) $(PY_DIR)/n7_gaffke_certificate.py --check-alpha 0.01

n7-certificate-005-check:
	$(PYTHON) $(PY_DIR)/n7_gaffke_certificate.py --check-alpha 0.05

n7-certificate-010-check:
	$(PYTHON) $(PY_DIR)/n7_gaffke_certificate.py --check-alpha 0.10

certificate-summary-check:
	@set -e; tmp=$$(mktemp); trap 'rm -f "$$tmp"' EXIT; \
	$(PYTHON) $(PY_DIR)/summarize_certificates.py --out $$tmp >/dev/null; \
	cmp $(CERT_DIR)/certificate-summary.json $$tmp

claim-manifest-check:
	$(PYTHON) $(PY_DIR)/validate_claim_manifest.py

paper:
	cd $(PAPER_DIR) && tectonic -X compile stringer.tex --keep-logs --keep-intermediates
	@! grep -Eq "Warning|Overfull|Underfull|undefined|multiply defined" $(PAPER_DIR)/stringer.log
	cd $(PAPER_DIR) && tectonic -X compile systematic-pps.tex --keep-logs --keep-intermediates
	@! grep -Eq "Warning|Overfull|Underfull|undefined|multiply defined" $(PAPER_DIR)/systematic-pps.log

reproduce: sync test public-corpus-check proof-check all-n-reduction-check tetrahedral-vertex-barrier-check five-coordinate-mixed-section-check five-coordinate-block-faces-check six-coordinate-repeated-top-check six-coordinate-unique-top-check six-coordinate-unique-top-l2-l3-check n4-weight-monotonicity-check n5-weight-monotonicity-check terminal-edge-bifurcation-check three-exponential-interior-check srswor-conditioning-check systematic-pps-check systematic-minimax-check systematic-disjoint-dp-check dirichlet-poissonization-check poisson-band-certificate-check poisson-band-calibration-check audit-benchmark-check simplex-cap-interface-check one-cap-all-n-check one-cap-certificate-check two-cap-certificate-check n3-formula-check n3-certificate-check n4-structure-check n4-certificate-check n5-structure-check n5-certificate-check n6-structure-check n6-certificate-check n7-certificate-check certificate-summary-check claim-manifest-check paper
	@echo "Core proofs, all-n reductions, complete four-, five-, and six-coordinate section-centroid barriers, terminal-edge bifurcation, three-exponential interior, SRSWOR conditioning, systematic-PPS failure and safeguard, finite-frame LP and disjoint-phase dynamic-program certificates, Dirichlet-Poissonization research certificates, Poisson simultaneous-band, scalar-calibration-path, decision-benchmark, simplex-cap challenge/solution interface, and uniform-multiplier certificates, analytic one-cap and certified two-cap checks, n=3 through n=7 Bernstein certificates, counterexample summaries, claim links, tests, and manuscript build passed."
