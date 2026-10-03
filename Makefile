.PHONY: install test pilot report

install:
	pip install -e "./python[dev]"

test:
	pytest python/ns_exploration/tests -v

pilot:
	python -m ns_exploration.experiments.pilot_campaign

sprint2-partial:
	python -m ns_exploration.experiments.sprint02_partial

sprint2-continued:
	python -m ns_exploration.experiments.sprint02_continued

sprint2-adjoint:
	python -m ns_exploration.experiments.sprint02_adjoint_ascent

sprint2-stress:
	python -m ns_exploration.experiments.sprint02_stress_c0001

sprint2-l0002:
	python -m ns_exploration.experiments.sprint02_l0002_cmaes

sprint2-l0003:
	python -m ns_exploration.experiments.sprint02_l0003_polish

sprint2-stress-c0002:
	python -m ns_exploration.experiments.sprint02_stress_c0002

sprint2-cs0001:
	python -m ns_exploration.experiments.sprint02_cs0001_l0005

sprint2-cs0001-attack:
	python -m ns_exploration.experiments.sprint02_cs0001_attack

sprint2-cert-stress:
	python -m ns_exploration.experiments.sprint02_cert_stress

sprint2-l0003-cert:
	python -m ns_exploration.experiments.sprint02_l0003_cert

sprint2-full-l0003-lean:
	python -m ns_exploration.experiments.sprint02_full_l0003_lean

sprint2-n16-rat:
	python -m ns_exploration.experiments.sprint02_n16_rat

sprint2-l0007:
	python -m ns_exploration.experiments.sprint02_l0007

sprint2-l0008:
	python -m ns_exploration.experiments.sprint02_l0008

sprint2-l0009:
	python -m ns_exploration.experiments.sprint02_l0009

sprint2-l0010:
	python -m ns_exploration.experiments.sprint02_l0010

sprint2-l0011:
	python -m ns_exploration.experiments.sprint02_l0011

sprint2-l0012:
	python -m ns_exploration.experiments.sprint02_l0012

sprint2-l0013:
	python -m ns_exploration.experiments.sprint02_l0013

sprint2-l0014:
	python -m ns_exploration.experiments.sprint02_l0014

sprint2-l0015:
	python -m ns_exploration.experiments.sprint02_l0015

sprint2-l0016:
	python -m ns_exploration.experiments.sprint02_l0016

sprint2-l0017:
	python -m ns_exploration.experiments.sprint02_l0017

sprint2-l0018:
	python -m ns_exploration.experiments.sprint02_l0018

sprint2-refute-c0002:
	python -m ns_exploration.experiments.sprint02_refute_c0002

sprint2-c0004:
	python -m ns_exploration.experiments.sprint02_c0004_n24

sprint2-l0020:
	python -m ns_exploration.experiments.sprint02_l0020

sprint2-l0021:
	python -m ns_exploration.experiments.sprint02_l0021

sprint2-l0022:
	python -m ns_exploration.experiments.sprint02_l0022

sprint2-l0023:
	python -m ns_exploration.experiments.sprint02_l0023

sprint2-l0024:
	python -m ns_exploration.experiments.sprint02_l0024

sprint2-l0025:
	python -m ns_exploration.experiments.sprint02_l0025

sprint2-l0026:
	python -m ns_exploration.experiments.sprint02_l0026

sprint2-l0027:
	python -m ns_exploration.experiments.sprint02_l0027

sprint2-l0028:
	python -m ns_exploration.experiments.sprint02_l0028

sprint2-l0029:
	python -m ns_exploration.experiments.sprint02_l0029

sprint2-l0030:
	python -m ns_exploration.experiments.sprint02_l0030

sprint2-l0031-0033:
	python -m ns_exploration.experiments.sprint02_l0031_0033

sprint2-l0034:
	python -m ns_exploration.experiments.sprint02_l0034

sprint2-l0035:
	python -m ns_exploration.experiments.sprint02_l0035

sprint2-l0036:
	python -m ns_exploration.experiments.sprint02_l0036

sprint2-l0037:
	python -m ns_exploration.experiments.sprint02_l0037

sprint2-l0038:
	python -m ns_exploration.experiments.sprint02_l0038

sprint2-l0039:
	python -m ns_exploration.experiments.sprint02_l0039

sprint2-l0040:
	python -m ns_exploration.experiments.sprint02_l0040

sprint2-l0041:
	python -m ns_exploration.experiments.sprint02_l0041

sprint2-l0042:
	python -m ns_exploration.experiments.sprint02_l0042

sprint2-l0043:
	python -m ns_exploration.experiments.sprint02_l0043

sprint2-l0044:
	python -m ns_exploration.experiments.sprint02_l0044

sprint2-l0045:
	python -m ns_exploration.experiments.sprint02_l0045

sprint2-l0046:
	python -m ns_exploration.experiments.sprint02_l0046

lean-build:
	cd lean/NSGalerkin && lake build

report:
	@echo "See reports/SPRINT_*.md"

docker-build:
	docker compose build

interval-test:
	julia --project=julia/ValidatedNS -e "using Pkg; Pkg.instantiate(); include(\"julia/ValidatedNS/test/runtests.jl\")"
