# AI-assisted exploration experiment

## Research question

What can be **built and checked** with AI-assisted software development when approaching a hard mathematical problem (3D Navier–Stokes regularity), and what failure modes appear when claims, certificates, and numerics are mixed without enough independent scrutiny?

## Personal / exploratory context

This repository is a personal lab (historical names: NS-MRL, bastardus2), not an institutional claim. It grew from ambitious attempts to push finite-dimensional bounds and tooling with heavy AI assistance. The public face is a **learning and integrity** case study, not a prize submission.

## Where AI assistance is evidenced

The archive contains large volumes of AI-shaped code, docs, reports, and tests. Concrete model names, versions, and prompt logs are **not systematically recorded** in-tree. Therefore:

| Claim | Stance |
|-------|--------|
| AI was used for coding / docs / exploration | Reasonable from project history and structure |
| Task X was done by model Y on date Z | **Not documented** — do not invent |

## What the repository lets you observe

- Working exploratory numerics and a small demo.
- Registries of conjectures and certificates with software verifiers.
- Documented corrections after invalid inequalities or verifier bugs (see `LESSONS_LEARNED.md`).

## What it does **not** measure

- Superiority of any AI model over another.
- Autonomous research capability in general.
- Progress toward solving the continuum Millennium problem (finite truncations are a different object).

## Runtime AI vs development AI

**There is no integration that calls an LLM or external model API during simulation or certificate verification.**  
“AI Lab” refers to **development assistance**, not an online AI solver. Do not add a runtime model dependency merely to justify the name.

## Human supervision

Degree of human review for each historical commit is **not documented** file-by-file. Treat unsupervised AI output as requiring verification. Failures in the code should **not** be automatically labeled “AI hallucination” unless provenance is known.

## No controlled evaluation

This repo alone is not an A/B study. No blind baselines, fixed budgets, or independent scoring panels are included.

## Optional protocol for future experiments (template)

Fill only with real data when you run a controlled trial:

| Field | Value |
|-------|--------|
| Task definition | |
| Model and version | |
| Prompt / context budget | |
| Allowed tools | |
| Success criterion | |
| Independent verifier | |
| Human review outcome | |
| Date / environment | |

Leave blank rather than inventing results.
