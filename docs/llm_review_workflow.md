# LLM Review Workflow

## Goal

Move from rule-based screening to a structured H&P-centered LLM adjudication workflow for possible immune checkpoint inhibitor toxicity.

## Recommended review strategy

Use H&P notes as the primary source for LLM review.

Why:

- H&P notes contain the strongest exposure and diagnostic narrative signal.
- They are less noisy than radiology for first-pass irAE adjudication.
- They usually contain the differential diagnosis, treatment context, and clinical reasoning needed to distinguish irAE from infection or progression.

## Script

Run:

```r
source("code/05_llm_review_prep.R")
```

## Inputs

The script expects these files to already exist:

- `output/note_extraction/encounter_level_rule_labels.csv`
- `output/note_extraction/all_notes_high_moderate_encounters.csv`
- optionally `output/checkpoint_irae_icu_epi/analysis_dataset.csv`

These are produced by:

```r
source("code/04_note_extraction.R")
```

## Outputs

Files are written to:

- `output/note_extraction/llm_review/`

Main outputs:

- `handp_top_notes_for_llm_review.csv`
  - top 3 H&P notes per flagged encounter
- `handp_llm_encounter_packets.csv`
  - encounter-level packet with concatenated note text and structured evidence columns
- `handp_llm_prompts.jsonl`
  - one prompt per encounter for local LLM inference
- `handp_manual_review_template.csv`
  - reviewer-friendly adjudication sheet
- `llm_review_prep_summary.csv`
  - summary counts for the prepared review set

## LLM task

For each encounter packet, the model should classify:

- whether ICI exposure is present
- whether exposure is recent or active
- whether irAE is suspected
- likely organ system
- strongest competing diagnosis
- whether steroids or rescue immunosuppression were given for suspected irAE
- overall encounter label:
  - `likely_ici_irae`
  - `possible_ici_irae`
  - `unlikely_ici_irae`

## Recommended human review loop

1. Run `code/04_note_extraction.R`.
2. Run `code/05_llm_review_prep.R`.
3. Run a local LLM over `handp_llm_prompts.jsonl`.
4. Merge model outputs back to `TEMP_ID` or `HAR`.
5. Review a validation sample in `handp_manual_review_template.csv`.
6. Refine the dictionary and prompt.
7. Lock a chart-reviewed training set for downstream modeling.

## Suggested first adjudication set

Start by manually reviewing:

- all encounters with both ICI evidence and irAE evidence
- all `high_priority_ici_irae_review` encounters
- a random sample of `possible_ici_irae` encounters
- a random sample of `unlikely_ici_irae` encounters for specificity checks

## Practical note

This stage is intended to create a high-quality labeled set for adjudication and eventual model development. It is not yet training a model directly.
