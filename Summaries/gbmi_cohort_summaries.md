GBMI Longitudinal Working Group Cohort Summaries
================

- [1 Set Up](#1-set-up)
- [2 Overview of scripts](#2-overview-of-scripts)
- [3 Overall Cohort Summary](#3-overall-cohort-summary)
- [4 Onset Summary](#4-onset-summary)
- [5 Progression Summary](#5-progression-summary)
- [6 Running the pipeline on all
  phenotypes](#6-running-the-pipeline-on-all-phenotypes)
- [7 Prepare files for sharing](#7-prepare-files-for-sharing)

<style>
.hl-line {
  background-color: #fff2a8;
  display: inline;
  box-shadow: 100vw 0 #fff2a8, -100vw 0 #fff2a8;
  clip-path: inset(0 -100vw);
}
</style>

<style>
pre.sourceCode code {
  white-space: pre-wrap;
  overflow-wrap: anywhere;
}
</style>

<script>
document.addEventListener("DOMContentLoaded", function() {
&#10;  document.querySelectorAll('[class*="highlight-lines-"]').forEach(function(el) {
&#10;    // Get specification from class, e.g. highlight-lines-8_9_12
    let cls = Array.from(el.classList)
      .find(x => x.startsWith("highlight-lines-"));
&#10;    if (!cls) return;
&#10;    let spec = cls.replace("highlight-lines-", "");
&#10;    // "_" separates line numbers
    let linesToHighlight = new Set(
      spec.split("_").map(Number)
    );
&#10;    // Find code element
    let code = el.tagName.toLowerCase() === "code"
      ? el
      : el.querySelector("code");
&#10;    if (!code) return;
&#10;    let lines = code.innerHTML.split("\n");
&#10;    code.innerHTML = lines.map(function(line, i) {
      if (linesToHighlight.has(i + 1)) {
        return '<span class="hl-line">' + line + '</span>';
      }
      return line;
    }).join("\n");
&#10;  });
&#10;});
</script>

# 1 Set Up

These scripts will compute descriptive summary statistics of the overall
biobank cohort, as well as each onset and progression phenotype.

The code is designed to use similar inputs as those provided for
separate scripts for PRS validation. Please ensure the following R
packages are available and installed.

```
install.packages(c("optparse", "data.table", "dplyr", "tidyr", "purrr", "survival", "broom"))
```

# 2 Overview of scripts

Every function call, regardless of `--mode`, needs three main flags:

- `--mode`: one of `cohort`, `onset`, or `progression`
- `--pheno-files`: one file, or a comma-separated list of files (see
  “Ancestry-specific vs. combined datasets” examples below). In
  multi-ancestry cohorts, summary statistics across all ancestry groups
  will also be generated. The script operates the same for biobanks and
  phenotypes with a single ancestry group.
- `--prefix`: where to write output, e.g. `results/mytrait/mytrait`. Any
  missing directories in the path are created automatically.

Everything else has a GBMI-convention default (`--id-col=IID`,
`--sex-col=sex`, `--event-col=event`, `--birthyear-col=birthyear`,
etc.), but should be overridden with `--*-col` flags to match your
actual column names. The script does **not** guess column names — a
mismatch raises a “column not found” error naming exactly which column
it looked for and couldn’t find.

## 2.1 Age vs. date columns

This applies to enrollment and follow-up in `cohort` mode. For each of
those two fields, the underlying age can be supplied two ways, and you
should pick **one** depending on what is most readily available:

- **Exact age** (e.g. `age_at_recruitment`): point the plain `-col` flag
  (`--enrollment-col` or `--followup-col`/`--followup-col-in-file`)
  directly at it. Used as-is, no conversion.
- **Date or year** (e.g. `year_at_recruitment`, or a `YYYY-MM-DD`
  string): point the matching `-date-col` flag (`--enrollment-date-col`
  or `--followup-date-col`/`--followup-date-col-in-file`) at it instead.
  This is converted to age using `--birthyear-col`. Bare 4-digit years
  are used directly; anything else is parsed as `YYYY-MM-DD`.

Don’t supply both flags for the same field. If enrollment or follow-up
live in a *separate* file from `--pheno-files` (via
`--enrollment-file`/`--followup-file`), the `-col`/`-date-col` flags
describe the column **in that separate file**, not in `--pheno-files` —
see the “Ancestry-specific datasets to be joined” example below.

`onset` and `progression` mode only use age since `--eventage-col` and
`--onsetage-col` are read straight from `--pheno-files` (used for
time-to-event GWAS), which should already have this information. For
progression specifically, if `--eventage-col` is not found in the
phenotype file, the script automatically derives it from
`--onsetage-col` + `--time-to-progression-col` (elapsed time since first
diagnosis) — see the progression examples below.

## 2.2 Ancestry-specific vs. combined datasets

`--pheno-files` accepts either:

- **One file** that already has every ancestry stacked together, with an
  ancestry column matching `--ancestry-col` — nothing further needed.
- **A comma-separated list** of per-ancestry files (no ancestry column
  required in each), paired with `--ancestry-labels` — a comma-separated
  list of labels, same order and same length as `--pheno-files` — which
  gets stamped onto each file before they’re combined.

If you’re looping over many phenotypes and not every one has a file for
every ancestry, add `--skip-missing-ancestries`: any file in the list
that doesn’t exist (and its paired label) is dropped with a warning,
rather than failing the whole run.

## 2.3 Histogram and KM masking for data privacy

Two independent flags control what gets suppressed before results leave
the biobank:

- **`--min-cell-count`** (default `0`, i.e. off): any histogram bin
  (age, birthyear, enrollment, follow-up, …) with a count below this
  threshold is set to `NA`, and a `minimum_cell_count_pass` column is
  added to that file so recipients can see exactly which bins were
  masked.
- **`--km-min`** (onset/progression only, default `0`, i.e. off):
  re-groups the Kaplan-Meier time axis so that every *published* time
  point has at least this many people at risk, experiencing the event,
  and censored — merging adjacent time bins together as needed. This
  protects the survival curve itself, which `--min-cell-count` does not
  touch (that flag only applies to the histogram tables).

Both are optional and independent — use one, both, or neither depending
on your biobank’s disclosure policy. The examples below marked “minimum
sample size” show both in action, including what the masked output looks
like compared to the unmasked version.

# 3 Overall Cohort Summary

This step only needs to be run one time. The key input information is

- sex
- ancestry
- birthyear
- age or date of enrollment
- age or date of last follow up

The input dataset can either be (1) combined file with all individuals
used for GWAS analysis or (2) list of ancestry-specific files. If
enrollment and/or last-followup information are contained in separate
files (each assumed to include all individuals), that can be combined in
this step.

## 3.1 Outputs

Every run in `cohort` mode writes the following files (all
gzip-compressed, tab-separated), stratified by ancestry x sex, with
`"ALL"` rollup rows across ancestry and/or sex:

- **`*_sample_sizes.txt.gz`** — `n` per ancestry and sex
- **`*_summaries.txt.gz`** — mean/median/sd of `birthyear`,
  `enrollment_age`, `enrollment_year`, `followup_age`, and
  `followup_time` (follow-up age minus enrollment age), per ancestry and
  sex
- **`*_birthyear_histogram.txt.gz`** — counts by ancestry, sex, birth
  year
- **`*_enrollmentAge_histogram.txt.gz`** — counts by ancestry, sex,
  enrollment age
- **`*_enrollmentYear_histogram.txt.gz`** — counts by ancestry, sex,
  enrollment year (`birthyear + enrollment_age`)
- **`*_followup_histogram.txt.gz`** — counts by ancestry, sex, total
  follow-up time (in years)

When `--min-cell-count` is set, each of the five histogram files above
gains a `minimum_cell_count_pass` column, and `n` is set to `NA`
wherever a bin fell below the threshold (see the “minimum sample size”
example and the masked-vs-unmasked comparison below).

Below are example scripts for different file availability scenarios
(example with UK Biobank)

## 3.2 Full dataset ready (age at enrollment / recruitment)

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=cohort \
--pheno-files=UKBB_cohort.txt \
--id-col=s \
--sex-col=sex \
--ancestry-col=pop \
--birthyear-col=birthyear \
--enrollment-col=age_at_recruitment \
--followup-col-in-file=age_death_or_lastvisit \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_cohort_summaries/UKBB_cohort
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## creating output directory  TTE_Summary/UKBB_cohort_summaries 
    ## RUN COMPLETE: yes

## 3.3 Full dataset ready (year of enrollment / recruitment)

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=cohort \
--pheno-files=UKBB_cohort.txt \
--id-col=s \
--sex-col=sex \
--ancestry-col=pop \
--birthyear-col=birthyear \
--enrollment-date-col=year_at_recruitment \
--followup-date-col-in-file=year_death_or_lastvisit \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_cohort_summaries2/UKBB_cohort
```

## 3.4 Full dataset ready (example with minimum sample size)

When `--min-cell-count` is specified for data privacy reasons (default
is 0), histogram bins below that minimum number are masked as NA. This
creates an additional column `minimum_cell_count_pass` in the histogram
files.

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=cohort \
--pheno-files=UKBB_cohort.txt \
--id-col=s \
--sex-col=sex \
--ancestry-col=pop \
--birthyear-col=birthyear \
--enrollment-col=age_at_recruitment \
--followup-col-in-file=age_death_or_lastvisit \
--male-value=1 --female-value=0 \
--min-cell-count=10 \
--prefix=TTE_Summary/UKBB_cohort_summaries_MIN/UKBB_cohort
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## creating output directory  TTE_Summary/UKBB_cohort_summaries_MIN 
    ## RUN COMPLETE: yes

## 3.5 Ancestry-specific datasets to be joined

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=cohort \
--pheno-files=UKBB_cohort_AFR.txt,UKBB_cohort_CSA.txt,UKBB_cohort_EUR.txt \
--ancestry-labels=AFR,CSA,EUR \
--id-col=s \
--sex-col=sex \
--ancestry-col=pop \
--birthyear-col=birthyear \
--enrollment-file=UKBB_enrollment_YEAR.txt --enrollment-date-col=year_at_recruitment --enrollment-id-col=s \
--followup-file=UKBB_followup.txt --followup-col-in-file=age_death_or_lastvisit --followup-id-col=s \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_cohort_summaries_ANC/UKBB_cohort
```

## 3.6 Example Results

![](gbmi_cohort_summaries_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# 4 Onset Summary

This step can be conducted with the phenotype files used to conduct TTE
GWAS. The key input information is

- sex
- ancestry
- birthyear
- `event`: onset case/control status
- `eventage`: age or date of diagnosis (cases) or last follow up
  (controls)

The key outputs are

- Summary statistics and histograms of age of diagnosis
- KM curves of disease onset

## 4.1 Outputs

Every run in `onset` mode writes (gzip-compressed, tab-separated,
overall and stratified by ancestry + sex, with `event` as an additional
case/control stratum):

- **`*_sample_sizes.txt.gz`** — `n` per ancestry, sex, `event` group
- **`*_summaries.txt.gz`** — mean/median/sd of age at event/censoring
  (`age`) and calendar year of event/censoring (`eventyear`), by
  ancestry, sex, `event`
- **`*_age_histogram.txt.gz`** — counts by ancestry, sex, age, `event`
- **`*_year_histogram.txt.gz`** — counts by ancestry, sex, calendar
  year, `event`
- **`*_km.txt.gz`** — Kaplan-Meier estimates of age-at-onset, one curve
  per ancestry and sex stratum: `time`, `n.risk`, `n.event`, `n.censor`,
  `estimate` (event-free probability), `std.error`,
  `conf.low`/`conf.high`, and `n_total`/`n_events_total` for that
  stratum

Note that `--min-cell-count` masks small counts (and adds
`minimum_cell_count_pass`) in the two histogram files above, and
`--km-min` re-bins `*_km.txt.gz` so every published time point meets the
minimum at-risk/event/censor counts — see the comparison plot at the end
of this section.

## 4.2 Combined-ancestry phenotype file

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=onset \
--pheno-files=UKBB_HYPTENSESS_pheno.txt \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=event \
--eventage-col=age \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESS/UKBB_HYPTENSESS 
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## creating output directory  TTE_Summary/UKBB_HYPTENSESS 
    ## RUN COMPLETE: yes

## 4.3 Combined-ancestry phenotype file (time unit conversion)

If event age is measured in days, please add the `--units` flag.
Allowable units are years (default), months, and days.

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=onset \
--pheno-files=UKBB_HYPTENSESS_pheno_days.txt \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=event \
--eventage-col=age \
--units='days' \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESS_DAY/UKBB_HYPTENSESS 
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## RUN COMPLETE: yes

## 4.4 Combined-ancestry phenotype file (with minimum sample sizes)

When `--km-min` is specified for data privacy reasons (default is 0),
the data table used to create KM curves is re-grouped so that time bins
include a sufficient number of samples.

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=onset \
--pheno-files=UKBB_HYPTENSESS_pheno.txt \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=event \
--eventage-col=age \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESS_MIN/UKBB_HYPTENSESS \
--km-min=10 \
--min-cell-count=10
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## creating output directory  TTE_Summary/UKBB_HYPTENSESS_MIN 
    ## RUN COMPLETE: yes

## 4.5 Ancestry-specific phenotype files

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=onset \
--pheno-files=UKBB_HYPTENSESS_pheno_AFR.txt,UKBB_HYPTENSESS_pheno_CSA.txt,UKBB_HYPTENSESS_pheno_EUR.txt \
--ancestry-labels=AFR,CSA,EUR \
--id-col=s \
--sex-col=sex \
--ancestry-col=pop \
--event-col=event \
--eventage-col=age \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESS_ANC/UKBB_HYPTENSESS
```

## 4.6 Ancestry-specific phenotype files (more general)

Feel free to keep a fixed list of ancestries if looping through multiple
phenotype files, as the script will skip over any missing files with the
flag `--skip-missing-ancestries`

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=onset \
--pheno-files=UKBB_HYPTENSESS_pheno_AFR.txt,UKBB_HYPTENSESS_pheno_CSA.txt,\
UKBB_HYPTENSESS_pheno_EUR.txt,UKBB_HYPTENSESS_pheno_OTHER.txt \
--ancestry-labels=AFR,CSA,EUR,OTHER \
--skip-missing-ancestries \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=event \
--eventage-col=age \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESS_ANC2/UKBB_HYPTENSESS
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## RUN COMPLETE: yes
    ## Warning message:
    ## In combine_ancestry_files(pheno_files, ancestry_col = ancestry_col,  :
    ##   1 of 4 ancestry file(s) not found and skipped: UKBB_HYPTENSESS_pheno_OTHER.txt (OTHER)

## 4.7 Example Results

![](gbmi_cohort_summaries_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](gbmi_cohort_summaries_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->![](gbmi_cohort_summaries_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

# 5 Progression Summary

Similar to onset summaries, this step can be conducted with the
phenotype files used to conduct TTE GWAS. The key input information is

- sex
- ancestry
- birthyear
- `event`: progression case/control status (may have a different name
  than onset)
- `eventage`: age or date of second diagnosis (progression cases) or
  last follow up (controls)
- `time-to-progression`: time from first diagnosis to second diagnosis
  or last follow up **optional, only to be used if event age not
  available**
- `onsetage`: age or date of first diagnosis

The key outputs are

- Summary statistics and histograms of baseline and progression age
- KM curves of disease progression by time from first diagnosis or
  absolute age

## 5.1 Outputs

Every run in `progression` mode writes (gzip-compressed, tab-separated,
stratified by ancestry, sex, `event`):

- **`*_sample_sizes.txt.gz`** — `n` per ancestry, sex, `event`
- **`*_summaries.txt.gz`** — mean/median/sd of onset age, event age, and
  event year, by ancestry, sex, `event`
- **`*_baselineAge_histogram.txt.gz`** — counts by ancestry, sex,
  age-at-first-diagnosis, `event`
- **`*_eventAge_histogram.txt.gz`** — counts by ancestry, sex,
  age-at-second-event/censoring, `event`
- **`*_[year/month]_histogram.txt.gz`** — counts by ancestry, sex,
  time-since-first-diagnosis, `event`
- **`*_age_km.txt.gz`** — Kaplan-Meier curve of progression on an
  absolute-age time scale
- **`*_[year/month]_km.txt.gz`** — Kaplan-Meier curve of progression on
  a time-since-first-diagnosis time axis

As above, `--min-cell-count` masks small counts and adds
`minimum_cell_count_pass` to the histogram files; `--km-min` re-bins the
`*_km.txt.gz` files.

## 5.2 Combined-ancestry phenotype file (age of progression column)

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=progression \
--pheno-files=UKBB_HYPTENSESStoCKD_pheno.txt \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=secondEvent \
--eventage-col=pheno2Age \
--onsetage-col=diagAge \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESStoCKD/UKBB_HYPTENSESStoCKD 
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## creating output directory  TTE_Summary/UKBB_HYPTENSESStoCKD 
    ## RUN COMPLETE: yes

## 5.3 Combined-ancestry phenotype file (time to progression column)

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=progression \
--pheno-files=UKBB_HYPTENSESStoCKD_pheno.txt \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=secondEvent \
--time-to-progression-col=secondTime \
--onsetage-col=diagAge \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESStoCKD2/UKBB_HYPTENSESStoCKD 
```

## 5.4 Combined-ancestry phenotype file (time unit conversion)

Time units are assumed to be in years for both endpoints, and thus two
units must be provided for `--units`. The default is `"years,years"`,
but below is an example with `onsetAge` measured in years and
`secondTime` (time from initial onset) measured in days.

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=progression \
--pheno-files=UKBB_HYPTENSESStoCKD_pheno_days.txt \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=secondEvent \
--time-to-progression-col=secondTime \
--onsetage-col=diagAge \
--units='years,days' \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESStoCKD_DAY/UKBB_HYPTENSESStoCKD 
```

    ## Warning messages:
    ## 1: package ‘optparse’ was built under R version 4.5.2 
    ## 2: package ‘data.table’ was built under R version 4.5.2 
    ## 3: package ‘dplyr’ was built under R version 4.5.2 
    ## 4: package ‘tidyr’ was built under R version 4.5.2 
    ## 5: package ‘purrr’ was built under R version 4.5.2 
    ## 6: package ‘survival’ was built under R version 4.5.2 
    ## 7: package ‘broom’ was built under R version 4.5.2 
    ## `age` not found in the phenotype file(s); deriving age at second event as `diagAge` + `secondTime` instead.
    ## RUN COMPLETE: yes

## 5.5 Ancestry-specific phenotype files

``` bash
Rscript gbmi_cohort_summaries.R \
--mode=progression \
--pheno-files=UKBB_HYPTENSESStoCKD_pheno_AFR.txt,UKBB_HYPTENSESStoCKD_pheno_CSA.txt,\
UKBB_HYPTENSESStoCKD_pheno_EUR.txt,UKBB_HYPTENSESStoCKD_pheno_OTHER.txt \
--ancestry-labels=AFR,CSA,EUR,OTHER \
--skip-missing-ancestries \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=secondEvent \
--eventage-col=pheno2Age \
--onsetage-col=diagAge \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_HYPTENSESStoCKD_ANC/UKBB_HYPTENSESStoCKD 
```

## 5.6 Example Results

![](gbmi_cohort_summaries_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->![](gbmi_cohort_summaries_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

# 6 Running the pipeline on all phenotypes

Summaries for all phenotypes can be computed and organized using the
following code (examples in bash and R below):

**NOTE**: Some summary statistics were not included in the final
meta-analysis following our QC procedures, e.g., due to small sample
size or observed inflation. Please ensure that those phenotype files are
not included in the summary analysis.

## 6.1 Bash (explicit list of phenotypes)

``` bash
onset_phenos=(HYPTENSESS T1D T2D CAD)

for onset_pheno in "${onset_phenos[@]}";; do

Rscript gbmi_cohort_summaries.R \
--mode=onset \
--pheno-files=UKBB_${onset_pheno}_pheno_AFR.txt,UKBB_${onset_pheno}_pheno_CSA.txt,UKBB_${onset_pheno}_pheno_EUR.txt,UKBB_${onset_pheno}_pheno_OTHER.txt \
--ancestry-labels=AFR,CSA,EUR,OTHER \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=event \
--eventage-col=age \
--skip-missing-ancestries \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_${onset_pheno}/UKBB_${onset_pheno}

done
```

``` bash
prog_phenos=(HYPTENSESStoCKD HYPTENSESStoMI HYPTENSESStoHEARTFAIL)

for prog_pheno in "${prog_phenos[@]}";; do

Rscript gbmi_cohort_summaries.R \
--mode=progression \
--pheno-files=UKBB_${prog_pheno}_pheno_AFR.txt,UKBB_${prog_pheno}_pheno_CSA.txt,UKBB_${prog_pheno}_pheno_EUR.txt,UKBB_${prog_pheno}_pheno_OTHER.txt \
--ancestry-labels=AFR,CSA,EUR,OTHER \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=secondEvent \
--eventage-col=pheno2Age \
--onsetage-col=diagAge \
--skip-missing-ancestries \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_${prog_pheno}/UKBB_${prog_pheno}

done
```

## 6.2 Bash (text file of phenotypes)

``` bash
while read onset_pheno; do

Rscript gbmi_cohort_summaries.R \
--mode=onset \
--pheno-files=UKBB_${onset_pheno}_pheno_AFR.txt,UKBB_${onset_pheno}_pheno_CSA.txt,UKBB_${onset_pheno}_pheno_EUR.txt,UKBB_${onset_pheno}_pheno_OTHER.txt \
--ancestry-labels=AFR,CSA,EUR,OTHER \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=event \
--eventage-col=age \
--skip-missing-ancestries \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_${onset_pheno}/UKBB_${onset_pheno}

done < onset_phenos.txt
```

``` bash
while read prog_pheno; do

Rscript gbmi_cohort_summaries.R \
--mode=progression \
--pheno-files=UKBB_${prog_pheno}_pheno_AFR.txt,UKBB_${prog_pheno}_pheno_CSA.txt,UKBB_${prog_pheno}_pheno_EUR.txt,UKBB_${prog_pheno}_pheno_OTHER.txt \
--ancestry-labels=AFR,CSA,EUR,OTHER \
--id-col=s \
--sex-col=sex \
--birthyear-col=birthyear \
--ancestry-col=pop \
--event-col=secondEvent \
--eventage-col=pheno2Age \
--onsetage-col=diagAge \
--skip-missing-ancestries \
--male-value=1 --female-value=0 \
--prefix=TTE_Summary/UKBB_${prog_pheno}/UKBB_${prog_pheno}

done < prog_phenos.txt
```

## 6.3 R

``` r
onset_phenos = c('HYPTENSESS', 'T1D', 'T2D', 'CAD')

for (onset_pheno in onset_phenos) {
  system2(command='Rscript',
          args=c(
            'gbmi_cohort_summaries.R',
            '--mode=onset',
            paste0('--pheno-files=UKBB_', onset_pheno, '_pheno_AFR.txt,UKBB_', onset_pheno, '_pheno_CSA.txt,UKBB_', onset_pheno, '_pheno_EUR.txt,UKBB_', onset_pheno, '_pheno_OTHER.txt'),
            '--ancestry-labels=AFR,CSA,EUR,OTHER',
            '--id-col=s',
            '--sex-col=sex',
            '--birthyear-col=birthyear',
            '--ancestry-col=pop',
            '--event-col=event',
            '--eventage-col=age',
            '--skip-missing-ancestries',
            '--male-value=1', '--female-value=0',
            paste0('--prefix=TTE_Summary/UKBB_', onset_pheno, '/UKBB_', onset_pheno)))
}
```

``` r
prog_phenos=c('HYPTENSESStoCKD', 'HYPTENSESStoMI', 'HYPTENSESStoHEARTFAIL')

for (prog_pheno in prog_phenos) {
  system2(command='Rscript',
          args=c(
            'gbmi_cohort_summaries.R',
            '--mode=progression',
            paste0('--pheno-files=UKBB_', prog_pheno, '_pheno_AFR.txt,UKBB_', prog_pheno, '_pheno_CSA.txt,UKBB_', prog_pheno, '_pheno_EUR.txt,UKBB_', prog_pheno, '_pheno_OTHER.txt'),
            '--ancestry-labels=AFR,CSA,EUR,OTHER',
            '--id-col=s',
            '--sex-col=sex',
            '--birthyear-col=birthyear',
            '--ancestry-col=pop',
            '--event-col=secondEvent',
            '--eventage-col=pheno2Age',
            '--onsetage-col=diagAge',
            '--skip-missing-ancestries',
            '--male-value=1', '--female-value=0',
            paste0('--prefix=TTE_Summary/UKBB_', prog_pheno, '/UKBB_', prog_pheno)))
}
```

# 7 Prepare files for sharing

All individual summary files are gzip compressed, and should be
organized within a TTE_Summary directory (as in the example script). For
final sharing, please upload the entire TTE_Summary directory to the
GBMI Google Cloud Bucket (e.g.,
`gsutil cp TTE_Summary/* gs://gbmi-[biobank]/TTE_Summary/` or
`gcloud storage rsync TTE_Summary gs://gbmi-[biobank]/TTE_Summary/`).

Thank you very much for your contributions to this project, and please
do not hesitate to contact us with any questions!

- Tony Chen <chentony@broadinstitute.org>
- Wei Zhou <wzhou@broadinstitute.org>
- Ida Surakka <isurakka@umich.edu>
