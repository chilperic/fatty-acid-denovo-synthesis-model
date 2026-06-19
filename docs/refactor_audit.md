# Refactor audit

Two uploaded FADNS archives were byte-identical after extraction:

- `semi-mechanistic-model-of-fatty-acids-de-novo-synthesis-main.zip`
- `fatty-acid-de-novo-synthesis-main (1).zip`

This repository uses one canonical copy.

Original issues addressed:

- notebook-only execution,
- `get_ipython()` side effects in exported Python files,
- repeated model code across variants,
- generated outputs mixed with source code,
- office lock files,
- duplicated data PDFs,
- no testable package boundary.

The refactor keeps sanitized legacy variant functions importable under `src/fadns_model/legacy_variants/`, while the reusable project utilities live in `src/fadns_model/`.
