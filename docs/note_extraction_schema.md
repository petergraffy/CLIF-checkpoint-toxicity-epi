# Note Extraction Schema

## Goal

Build a hybrid rules-plus-LLM workflow to identify evidence of prior immune checkpoint inhibitor exposure, suspected immune-related adverse events (irAEs), competing diagnoses, treatment response, and radiographic support from H&P notes and radiology reports.

## Input files

### H&P notes

Expected columns from `H&P_VAR.csv`:

- `MRN`
- `HAR`
- `PAT_ENC_CSN_ID`
- `ENC_DATE`
- `CONTACT_DATE`
- `NOTE_ID`
- `NOTE_CSN_ID`
- `NOTE_TYPE`
- `NOTE_STATUS`
- `NOTE_ENC_DATE`
- `NOTE_DTTM`
- `SERVICE_DTTM`
- `AUTHOR_NAME`
- `AUTHOR_SERV`
- `NOTE_TEXT`

### Radiology reports

Expected columns from `IMG_VAR.csv`:

- `MRN`
- `HAR`
- `PROC_ID`
- `PROC_CODE`
- `PROC_NAME`
- `ORD_PROC_ID`
- `ORDER_ID`
- `ACCESSION_NUM`
- timing fields such as `ORDERING_DTTM`, `BEGIN_EXAM_DTTM`, `END_EXAM_DTTM`, `FINALIZING_DTTM`
- `NOTE_CSN_ID`
- `NOTE_ID`
- `NOTE_TEXT`

## Rule-based note-level extraction targets

### Exposure

- explicit ICI drug mention
- general immunotherapy or checkpoint inhibitor mention
- potential recency signal such as "last dose", "cycle", "received on"

### irAE suspicion

- general irAE terminology
- suspected immunotherapy toxicity
- organ-specific toxicity:
  - pneumonitis
  - myocarditis
  - hepatitis
  - colitis
  - nephritis
  - endocrine toxicity
  - neurologic toxicity
  - HLH-like toxicity

### Competing explanations

- infection or pneumonia
- cancer progression
- pulmonary edema or heart failure
- aspiration

### Treatment response

- corticosteroids
- second-line rescue immunosuppression

### Qualifiers

- favoring language
- uncertainty language

### Radiology-only support

- pneumonitis-compatible pattern
- infection-compatible pattern
- edema-compatible pattern
- progression-compatible pattern

## Encounter-level label set

Primary grouping variable:

- `HAR`

Recommended encounter-level labels:

- `encounter_any_ici_exposure`
- `encounter_explicit_ici_agent`
- `encounter_any_irae_signal`
- `encounter_pneumonitis_signal`
- `encounter_cardiac_signal`
- `encounter_hepatitis_signal`
- `encounter_colitis_signal`
- `encounter_renal_signal`
- `encounter_endocrine_signal`
- `encounter_neurologic_signal`
- `encounter_hlh_signal`
- `encounter_competing_infection`
- `encounter_competing_progression`
- `encounter_competing_edema`
- `encounter_competing_aspiration`
- `encounter_steroid_signal`
- `encounter_rescue_signal`
- `encounter_radiology_supportive`
- `encounter_llm_priority`

Rule-based encounter summary categories:

- `high_priority_ici_irae_review`
- `moderate_priority_possible_ici_irae`
- `competing_diagnosis_only`
- `low_signal`

## LLM review queue

The rule-based pass should output a smaller encounter-focused queue for local LLM review. Recommended fields:

- `MRN`
- `HAR`
- `NOTE_ID`
- `note_source`
- `note_datetime`
- `NOTE_TYPE` or `PROC_NAME`
- `rule_score`
- `llm_priority_note`
- `evidence_snippet_ici`
- `evidence_snippet_irae`
- `evidence_snippet_competing`
- `NOTE_TEXT`

## Recommended first LLM tasks

For each note:

- confirm whether prior ICI exposure is present
- determine whether irAE is explicitly suspected, favored, or ruled out
- classify likely organ system
- identify strongest competing diagnosis
- identify whether clinicians initiated steroids or rescue immunosuppression for suspected irAE
- return structured JSON plus a brief evidence rationale
