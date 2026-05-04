# CLIF Checkpoint Inhibitor–Associated Critical Illness Phenotyping

## CLIF VERSION

2.1.0

## Objective

This project identifies adult ICU hospitalizations among patients with cancer present on admission (POA) and phenotypes a subset with clinical features suggestive of severe immune checkpoint inhibitor–associated toxicity or immune-related adverse events (irAEs). Because checkpoint inhibitors are typically administered outside the ICU and may not be reliably represented in inpatient medication administration records, this project does **not** define exposure using checkpoint inhibitor medication tables. Instead, it uses a reverse-phenotyping approach within CLIF to identify ICU admissions with acute inflammatory, organ-specific, and treatment-response patterns compatible with severe irAE syndromes.

The primary goals are to:
1. Identify a broad adult cancer ICU cohort.
2. Characterize acute ICU syndromes consistent with severe irAE-like illness using biomarkers, vitals, organ support, and steroid/rescue treatment.
3. Classify patients into possible, probable, and high-confidence checkpoint inhibitor–associated ICU phenotypes for downstream validation and epidemiologic study.
4. Estimate whether the burden of suspected severe irAE ICU presentations is changing over time.
5. Evaluate whether county-level air pollution burden is associated with suspected irAE ICU phenotypes and worse outcomes.

## Required CLIF tables and fields

Please refer to the [CLIF data dictionary](https://clif-icu.com/data-dictionary), [CLIF Tools](https://clif-icu.com/tools), [ETL Guide](https://clif-icu.com/etl-guide), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for more information on constructing the required tables and fields.

The following tables are required:

1. **hospitalization**
   - `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `admission_type_name`, `admission_type_category`, `discharge_name`, `discharge_category`
   - Used to define hospitalization-level cohort membership and adult status.

2. **adt**
   - `hospitalization_id`, `in_dttm`, `out_dttm`, `location_category`
   - Used to identify ICU stays and derive the first ICU admission time (`first_icu_in`), last ICU exit time, and number of ICU segments.

3. **hospital_diagnosis**
   - `hospitalization_id`, `diagnosis_code`, `diagnosis_code_format`, `poa_present`
   - Used to identify malignancy present on admission and to derive supportive diagnosis-based organ flags. Hospital diagnosis is used as background/supportive context rather than the primary acute syndrome engine.

4. **labs**
   - `hospitalization_id`, `lab_result_dttm`, `lab_collect_dttm`, `lab_category`, `lab_name`, `lab_value`
   - Used to derive inflammatory, cardiac, hepatic, renal, hematologic, and hyperinflammatory features.
   - Relevant lab concepts include ferritin, CRP, ESR, troponin, AST, ALT, bilirubin, creatinine, lactate, WBC, eosinophils, and platelets.

5. **vitals**
   - `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value`
   - Used to derive acute physiologic severity around ICU entry.
   - Relevant `vital_category` values include `heart_rate`, `respiratory_rate`, `temp_c`, `map`, and `spo2`.

6. **respiratory_support**
   - `hospitalization_id`, `recorded_dttm`, `device_category`, `fio2_set`
   - Used to identify invasive mechanical ventilation, high-flow nasal cannula, noninvasive positive pressure ventilation, and oxygen requirement severity.

7. **crrt_therapy**
   - `hospitalization_id`, `recorded_dttm`
   - Used to identify severe renal failure requiring continuous renal replacement therapy.

8. **medication_admin_continuous**
   - `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`
   - Used to identify vasopressor exposure as a marker of shock severity.
   - Relevant medications include norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, dobutamine, milrinone, and angiotensin.

9. **medication_admin_intermittent**
   - `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `mar_action_group`
   - Used to identify steroids and rescue immunosuppression given shortly before or after ICU entry.
   - Relevant steroids include methylprednisolone, prednisone, prednisolone, dexamethasone, and hydrocortisone.
   - Relevant rescue immunosuppressants include mycophenolate, infliximab, tocilizumab, abatacept, rituximab, cyclophosphamide, IVIG, and immune globulin.

### Notes on project design

- **Checkpoint inhibitor medications are not used to define exposure**, because many relevant patients received immunotherapy before the ICU admission and those administrations may not appear in inpatient CLIF medication records.
- **Patient diagnosis is not required** for this workflow. Because `patient_diagnosis` may currently function as a concept table rather than a reliably populated encounter-level source, this project does not rely on it for acute syndrome timing.
- The phenotype is therefore driven primarily by:
  - cancer POA,
  - ICU timing,
  - biomarkers,
  - organ support,
  - steroids/rescue treatment,
  - and supportive diagnosis codes.

## Cohort identification

### Broad cohort definition

The broad study cohort includes:
1. Adult hospitalizations (`age_at_admission >= 18`)
2. With at least one ICU stay identified from `adt` where `location_category == "icu"`
3. With cancer present on admission (`poa_present = 1`) based on `hospital_diagnosis`

### Cancer definition

Cancer POA is defined using broad ICD-10-CM and ICD-9-CM malignant neoplasm code families in `hospital_diagnosis`, including:
- ICD-10 malignant neoplasm ranges broadly spanning `C00`–`C96`, including hematologic malignancies
- Personal/history-of-cancer style codes such as `Z85` or `V10` where available in the hospitalization diagnosis data

This broad oncology anchor is intended to capture ICU admissions among patients with an underlying malignancy, rather than restrict the cohort to any single cancer subtype.

### Exclusions

Current exclusions are minimal:
- age < 18
- no ICU segment
- no cancer diagnosis present on admission

Future versions may add exclusions for alternative ICU syndromes that would strongly compete with an irAE-like interpretation, such as trauma, perioperative admissions, or clearly noninflammatory critical illness phenotypes.

## Epidemiology aims

### Aim 1. Burden and severity over time

- Estimate annual rates of possible, probable, and high-confidence irAE-like ICU admissions among all cancer ICU admissions.
- Compare severity and outcomes for suspected irAE cases versus other cancer ICU admissions using:
  - hospital mortality
  - ICU and hospital length of stay
  - invasive mechanical ventilation
  - vasopressors
  - CRRT
  - discharge to hospice or facility

### Aim 2. Air pollution and social context

- Link CLIF hospitalization geography (`county_code`, `zipcode_five_digit`, or `census_tract` when available) to county-level PM2.5, NO2, SVI, and ACS covariates.
- Examine whether higher chronic pollution burden is associated with:
  - the suspected irAE phenotype among cancer ICU admissions
  - more severe organ support needs
  - mortality or prolonged ICU stay among suspected irAE cases

Because CLIF v2.1 includes encounter geography in `hospitalization`, this linkage can be done without storing full address data in the repo.

## TL;DR phenotype definitions

This project uses a **reverse phenotype** for checkpoint inhibitor–associated critical illness. Since prior outpatient immunotherapy exposure cannot be reliably observed in CLIF alone, the workflow identifies cancer ICU admissions with acute organ injury, inflammatory signatures, severe support needs, and steroid/rescue treatment patterns compatible with severe irAEs.

### Supportive diagnosis flags

Hospital diagnosis codes are used only as **supportive context**, not as the primary engine of the phenotype. Supportive diagnosis groups include:
- pneumonitis-like
- myocarditis/pericarditis-like
- hepatitis-like
- colitis-like
- nephritis/AKI-like
- endocrine toxicity-like
- neurologic inflammatory toxicity-like
- HLH-like

These are summarized into `organ_dx_count`.

### Severity features

Acute severity is derived near ICU entry and includes:
- **Respiratory support:** IMV, HFNC, NIPPV, high FiO2
- **Shock support:** vasopressors
- **Renal support:** CRRT
- **Vitals:** fever, low MAP, hypoxemia

These are summarized in `severe_support_count`.

### Biomarker features

Acute biologic signal is defined using labs obtained around ICU entry:
- **Inflammatory markers:** ferritin, CRP, ESR
- **Cardiac injury:** troponin
- **Hepatic injury:** AST, ALT, bilirubin
- **Renal injury:** creatinine
- **Shock/perfusion:** lactate
- **Hematologic/inflammatory support:** WBC, eosinophils, platelets

Derived summary flags include:
- `inflam_marker_positive`
- `hepatitis_lab_positive`
- `myocarditis_lab_positive`
- `renal_lab_positive`
- `cytokine_storm_like`

These are summarized in `immune_signal_count`.

### Steroid and rescue treatment features

Treatment-response features are defined using intermittent medication administration around ICU entry:
- `any_steroid`
- `any_rescue_immunosuppression`

These are important because many suspected severe irAE patients will receive high-dose steroids and, in some cases, second-line immunosuppressive therapy.

## Organ-specific irAE-like modules

The core phenotype is driven by five biomarker/treatment-based modules.

### 1. Pulmonary irAE-like
Flags ICU admissions with:
- severe respiratory support or hypoxemia
- and steroids

Operationally, this captures respiratory failure patterns compatible with pneumonitis-like toxicity.

### 2. Cardiac irAE-like
Flags ICU admissions with:
- elevated troponin
- and vasopressors and/or elevated lactate
- and steroids

Operationally, this captures myocarditis-like or inflammatory cardiogenic shock–like presentations.

### 3. Hepatic irAE-like
Flags ICU admissions with:
- marked transaminitis and/or hyperbilirubinemia
- and steroids

Operationally, this captures hepatitis-like toxicity.

### 4. Renal irAE-like
Flags ICU admissions with:
- elevated creatinine and/or CRRT
- and steroids

Operationally, this captures nephritis-like or severe renal inflammatory injury phenotypes.

### 5. Hyperinflammatory irAE-like
Flags ICU admissions with:
- elevated ferritin, CRP, ESR, or cytokine-storm-like lab pattern
- plus fever or severe organ support
- and steroids

Operationally, this captures HLH-like, cytokine-storm-like, or generalized hyperinflammatory immune toxicity patterns.

These are summarized in `n_irae_modules`.

## Final tiered phenotypes

### Possible checkpoint inhibitor irAE ICU
Defined as either:
- at least one biomarker/treatment-based irAE-like module, **or**
- supportive diagnosis evidence plus steroids plus severe support

This is the broadest phenotype and is intended for sensitivity.

### Probable checkpoint inhibitor irAE ICU
Defined as:
- at least one irAE-like module
- plus invasive mechanical ventilation, vasopressors, or CRRT

This reflects patients with a more severe and clinically persuasive phenotype.

### High-confidence checkpoint inhibitor irAE ICU
Defined as either:
- at least two distinct irAE-like modules plus major organ support, **or**
- rescue immunosuppression plus major organ support

This is the most specific phenotype and is intended to identify patients with the strongest computable evidence of severe immune-mediated toxicity.

## Expected Results

This project produces:
1. A broad adult cancer ICU cohort
2. A hospitalization-level feature table containing supportive diagnosis flags, biomarker features, vital sign features, organ support features, treatment-response features, and final phenotype classifications
3. A filtered cohort of suspected checkpoint inhibitor–associated critical illness hospitalizations
4. Summary counts for the broad cohort and each phenotype tier

Expected output files include:
- `output/checkpoint_irae_icu/01_cancer_icu_broad.csv`
- `output/checkpoint_irae_icu/02_cancer_icu_irae_features.csv`
- `output/checkpoint_irae_icu/03_cancer_icu_checkpoint_suspected.csv`
- `output/checkpoint_irae_icu/summary_counts.csv`

The final project results should be saved in the `output/final` directory if adapted into a production CLIF project structure.

## Detailed Instructions for Running the Project

### 1. Clone the repo and create `config/config.json`

Follow [config/README.md](config/README.md) and create `config/config.json`.

At minimum, most sites will need:

```json
{
  "site_name": "SITE",
  "clif_dir": "/path/to/clif-2.1.0",
  "output_dir": "/path/to/repo/output/checkpoint_irae_icu",
  "analysis_dir": "/path/to/repo/output/checkpoint_irae_icu_epi",
  "ollama_model": "qwen3:235b",
  "ollama_temperature": 0
}
```

### 2. Install the required R packages

Typical package requirements for this project are:

- `tidyverse`
- `lubridate`
- `arrow`
- `readr`
- `stringr`
- `broom`
- `jsonlite`

If you are using `renv`, restore the environment from your site-specific setup. If not, install packages directly in R:

```r
install.packages(c(
  "tidyverse",
  "lubridate",
  "arrow",
  "readr",
  "stringr",
  "broom",
  "jsonlite"
))
```

### 3. Install Ollama explicitly

The LLM review stage uses a local Ollama runtime.

#### macOS

Install Ollama with Homebrew:

```bash
brew install ollama
```

Or install the desktop app from the official site:

- [Ollama downloads](https://ollama.com/download)

After installation, verify that Ollama is available:

```bash
ollama --version
```

Start the Ollama service if needed:

```bash
ollama serve
```

In a second terminal, pull the default review model:

```bash
ollama pull qwen3:235b
```

If you prefer to use the local API rather than the CLI, the scripts will also work against the default Ollama host:

- `http://127.0.0.1:11434`

### 4. Place note files in the repo

For the note workflow, place exported note files here:

- `notes/HP_NOTES.csv`
- `notes/RAD_NOTES.csv`

The note extraction scripts assume those exact filenames unless you override them in `config/config.json`.

### 5. Run the pipeline in order

The pipeline is now designed to run stepwise from saved outputs.

```r
source("code/01_cohort.R")
source("code/02_analysis.R")
source("code/03_air_pollution_analysis.R")
source("code/04_note_extraction.R")
source("code/05_llm_review_prep.R")
source("code/06_ollama_llm_review.R")
```

### 6. What each step produces

`code/01_cohort.R`
- writes:
  - `output/checkpoint_irae_icu/01_cancer_icu_broad.csv`
  - `output/checkpoint_irae_icu/02_cancer_icu_irae_features.csv`
  - `output/checkpoint_irae_icu/03_cancer_icu_checkpoint_suspected.csv`
  - `output/checkpoint_irae_icu/summary_counts.csv`

`code/02_analysis.R`
- reads the saved cohort outputs from `01`
- writes:
  - `output/checkpoint_irae_icu_epi/analysis_dataset.csv`
  - denominator and QC tables

`code/03_air_pollution_analysis.R`
- reads `analysis_dataset.csv`
- links county-level pollution files
- writes analytical tables and figures

`code/04_note_extraction.R`
- reads `notes/HP_NOTES.csv` and `notes/RAD_NOTES.csv`
- writes note-level, encounter-level, and LLM queue files to `output/note_extraction/`

`code/05_llm_review_prep.R`
- builds H&P-centered encounter packets for LLM review
- writes prompt-ready files to `output/note_extraction/llm_review/`

`code/06_ollama_llm_review.R`
- runs local Ollama review over the H&P encounter packets
- writes:
  - `ollama_handp_raw_responses.csv`
  - `ollama_handp_parsed_outputs.csv`
  - `ollama_handp_review_merged.csv`
  - `ollama_handp_review_summary.csv`

### 7. Practical notes for federated sites

- Sites do not need to run the air pollution step to run the note extraction and LLM review steps.
- Sites do need `analysis_dataset.csv` for the note time-window alignment in `code/04_note_extraction.R`.
- If a site wants to use a different Ollama model, set `ollama_model` in `config/config.json`.
- The recommended default for one-shot extraction is `qwen3:235b` with `ollama_temperature = 0`.
