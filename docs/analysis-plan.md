# ICU Checkpoint Toxicity Epidemiology Analysis Plan

## Main framing

This project is set up around a reverse phenotype for severe immune checkpoint inhibitor (ICI) toxicity in ICU patients with cancer. That is a pragmatic design for CLIF because outpatient ICI administrations are often not captured in inpatient medication tables, while the downstream severe syndrome often is.

## Aim 1. Is severe suspected ICI toxicity becoming more common?

Primary denominator:
- Adult cancer ICU hospitalizations

Primary case definitions:
- `possible_checkpoint_irae_icu`
- `probable_checkpoint_irae_icu`
- `high_confidence_checkpoint_irae_icu`

Recommended burden metrics:
- Annual count
- Annual rate per 1,000 cancer ICU admissions
- Joinpoint-free temporal trend estimates using calendar year as a continuous term first

Recommended severity outcomes:
- `hospital_mortality`
- `any_imv`
- `any_vasopressor`
- `any_crrt`
- `prolonged_icu_los`
- `prolonged_hospital_los`
- hospice discharge

Key caution:
- A rising phenotype rate could reflect more true toxicity, higher ICI uptake, better recognition, or changing coding and treatment patterns. Treat this as severe suspected irAE burden unless outpatient oncology exposure data can be linked.

## Aim 2. Is air pollution associated with suspected ICI toxicity or worse outcomes?

Minimum linkage variable from CLIF v2.1:
- `county_code` from `hospitalization`

Current repo inputs already support:
- county-year PM2.5
- county-year NO2
- county-year SVI
- county/month pollution files for higher-resolution future work

Recommended first-pass analyses:
- Among all cancer ICU admissions: association of PM2.5 and NO2 with `probable_checkpoint_irae_icu`
- Among suspected irAE cases: association of PM2.5 and NO2 with mortality and major organ support

Minimum adjustment set for exploratory models:
- age
- calendar year
- broad cancer type
- SVI or ACS socioeconomic variables

Better adjustment set if available from CLIF:
- sex
- race/ethnicity
- transfer status
- code status
- baseline comorbidity burden

## Immediate next build steps

1. Save `analysis_dataset.csv` from `code/02_analysis.R`.
2. Confirm column names in the county-level PM2.5, NO2, and SVI files.
3. Run `code/03_air_pollution_analysis.R` and inspect linkage yield by county and year.
4. Add the `patient` table to the cohort pipeline so demographics are available for adjusted models.
5. If feasible, add an external oncology exposure source to separate "true ICI-exposed" cases from the broader reverse phenotype.
