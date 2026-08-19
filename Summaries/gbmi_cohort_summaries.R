#!/usr/bin/env Rscript
# ==============================================================================
# gbmi_cohort_summaries.R - Basic cohort / onset / progression summary stats
#                            for cross-biobank GWAS meta-analysis requests
# ==============================================================================
# Purpose:
#   Produce the two things requested from each biobank analyst:
#
#   1) overall_summary() - combined cohort demographics: ancestry, sex, ages of
#      enrollment and last follow-up. Handles the gaps biobanks typically hand
#      you:
#        - ancestry-specific phenotype files that need to be stacked into one
#          table with an `ancestry` column (combine_ancestry_files())
#        - enrollment age/date living in a separate file entirely
#          (join_enrollment())
#        - last-followup age/date living in a separate file entirely, which
#          may or may not be the same file as enrollment (join_followup())
#        - sex coded as a site-specific 0/1 (or other) convention instead of
#          "Male"/"Female" (recode_sex())
#
#   2) onset_summary() / prog_summary() - event-time summaries (sample sizes,
#      age/year histograms, KM curves) computed directly from the existing
#      ancestry-specific phenotype files, again stacked via
#      combine_ancestry_files() and sex-recoded via recode_sex().
#
# This intentionally does NOT reproduce the full 00_run_pipeline_gbmi.sh
# final/ + intermediate/ output scaffolding or its run-status/manifest
# machinery - it's meant to be a low-overhead companion script, not a second
# pipeline. What it *does* borrow, because skipping them would be a real
# disclosure/consistency risk:
#   - the "check package versions, don't silently install" policy used by
#     01-06 (biobank compute nodes often have no internet access, and silent
#     installs let sites drift onto different package versions)
#   - prs_risk_utils.R's suppress_small_cells(), sourced the same way 05/06
#     do it (via PRS_SCRIPT_DIR / --file), instead of a second, inconsistent
#     count-suppression rule living only in this script
#
# Usage (called directly with Rscript, or sourced interactively - see the
# runnable example guarded by is_main() at the bottom):
#
#   Rscript gbmi_cohort_summaries.R \
#     --mode=cohort \
#     --pheno-files=eur_demog.txt,afr_demog.txt,eas_demog.txt \
#     --ancestry-labels=EUR,AFR,EAS \
#     --birthyear-col=birthyear \
#     --enrollment-file=recruitment.txt --enrollment-col=age_at_recruitment \
#     --followup-file=followup.txt --followup-col-in-file=age_at_last_followup \
#     --male-value=1 --female-value=0 \
#     --prefix=mybiobank_cohort
#
#   Rscript gbmi_cohort_summaries.R \
#     --mode=onset \
#     --pheno-files=eur_pheno_PHECODE.txt,afr_pheno_PHECODE.txt \
#     --ancestry-labels=EUR,AFR \
#     --event-col=event \ 
#     --eventage-col=age \ 
#     --male-value=1 --female-value=0 \
#     --prefix=mybiobank_PHECODE --km-min=0
#
#   # --eventage-col recorded in days instead of years (e.g. some EHR-derived
#   # extracts) - convert with --units=days; also accepts --units=months.
#   # Default is --units=years, i.e. a no-op, so existing calls without this
#   # flag behave exactly as before.
#   Rscript gbmi_cohort_summaries.R \
#     --mode=onset \
#     --pheno-files=eur_pheno_PHECODE.txt,afr_pheno_PHECODE.txt \
#     --ancestry-labels=EUR,AFR \
#     --event-col=event \
#     --eventage-col=age_days --units=days \
#     --male-value=1 --female-value=0 \
#     --prefix=mybiobank_PHECODE --km-min=0
#
#   Rscript gbmi_cohort_summaries.R \
#     --mode=progression \
#     --pheno-files=eur_pheno_PHECODE1toPHECODE2.txt,afr_pheno_PHECODE1toPHECODE2.txt \
#     --ancestry-labels=EUR,AFR \
#     --event-col=secondEvent \
#     --eventage-col=pheno2Age \
#     --onsetage-col=diagAge \
#     --male-value=1 --female-value=0 \
#     --prefix=mybiobank_PHECODE1toPHECODE2 --km-min=0
#
#   # same, but the phenotype file only has time-to-progression (e.g.
#   # secondTime, in years), not an absolute age at second event - eventAge
#   # is derived automatically as onsetAge + secondTime when --eventage-col
#   # is missing from the file (see prog_summary())
#   Rscript gbmi_cohort_summaries.R \
#     --mode=progression \
#     --pheno-files=eur_pheno_PHECODE1toPHECODE2.txt,afr_pheno_PHECODE1toPHECODE2.txt \
#     --ancestry-labels=EUR,AFR \
#     --event-col=secondEvent \
#     --eventage-col=pheno2Age --time-to-progression-col=secondTime \
#     --onsetage-col=diagAge \
#     --male-value=1 --female-value=0 \
#     --prefix=mybiobank_PHECODE1toPHECODE2 --km-min=0
#
#   # progression mode, but the site's diagAge/pheno2Age (and secondTime,
#   # if used) are all recorded in months instead of years - one --units
#   # flag converts baseline (--onsetage-col), event (--eventage-col), and
#   # the time-to-progression fallback (--time-to-progression-col) together,
#   # since they all come from the same file in the same units
#   Rscript gbmi_cohort_summaries.R \
#     --mode=progression \
#     --pheno-files=eur_pheno_PHECODE1toPHECODE2.txt,afr_pheno_PHECODE1toPHECODE2.txt \
#     --ancestry-labels=EUR,AFR \
#     --event-col=secondEvent \
#     --eventage-col=pheno2Age_months --onsetage-col=diagAge_months --units=months \
#     --male-value=1 --female-value=0 \
#     --prefix=mybiobank_PHECODE1toPHECODE2 --km-min=0
#
# Author: (adapted from Ying Wang's Longitudinal-PRS pipeline conventions)
# ==============================================================================

# ==============================================================================
# PACKAGE MANAGEMENT (same policy as 01-06: check, don't silently install)
# ==============================================================================
load_package <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    if (identical(Sys.getenv("PRS_ALLOW_INSTALL"), "1")) {
      cat("Installing missing package (PRS_ALLOW_INSTALL=1):", pkg, "\n")
      install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
    } else {
      stop("Required package not installed: ", pkg,
           "\n  Install it, then re-run:\n    install.packages(\"", pkg, "\")",
           "\n  (or set PRS_ALLOW_INSTALL=1 to allow automatic installation)")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
required_packages <- c("optparse", "data.table", "dplyr", "tidyr", "purrr", "survival", "broom")
invisible(lapply(required_packages, load_package))

# ==============================================================================
# Source shared helper (same PRS_SCRIPT_DIR / --file resolution as 05/06,
# so this script and the main pipeline agree on what "disclosure-safe" means)
# ==============================================================================
resolve_script_dir <- function() {
  env_dir <- Sys.getenv("PRS_SCRIPT_DIR", unset = "")
  if (nzchar(env_dir)) return(normalizePath(env_dir, winslash = "/", mustWork = TRUE))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  script_path <- sub("^--file=", "", file_arg[[1]])
  script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
  dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
}
.this_dir <- resolve_script_dir()
.utils <- file.path(.this_dir, "prs_risk_utils.R")
if (file.exists(.utils)) {
  source(.utils)
} else {
  # Fall back so this script still runs standalone (e.g. copied to a site
  # that only wants the cohort summaries, not the full pipeline), but warn
  # loudly since the fallback below is NOT identical to the pipeline's rule.
  warning("prs_risk_utils.R not found next to this script; using a local, ",
          "less-tested suppress_small_cells(). Copy prs_risk_utils.R beside ",
          "this script to match the main pipeline's disclosure rule exactly.")
  suppress_small_cells <- function(x, count_cols, threshold = 10L, estimate_cols = character()) {
    x <- as.data.table(x)
    count_cols    <- intersect(count_cols, names(x))
    estimate_cols <- intersect(estimate_cols, names(x))
    if (!length(count_cols)) { x[, minimum_cell_count_pass := TRUE]; return(x[]) }
    small <- Reduce(`|`, lapply(count_cols, function(cc) {
      v <- suppressWarnings(as.numeric(x[[cc]])); is.finite(v) & v > 0 & v < threshold
    }))
    small[is.na(small)] <- FALSE
    x[, minimum_cell_count_pass := !small]
    if (any(small)) for (cc in unique(c(count_cols, estimate_cols))) x[small, (cc) := NA]
    x[]
  }
}

# ==============================================================================
# combine_ancestry_files() - stack ancestry-specific files into one table with
# a proper `ancestry` column.
#
# Biobanks are commonly asked to produce one phenotype/demographic file PER
# ANCESTRY GROUP. This handles both flavors analysts send back:
#   - files that already contain an ancestry column (nothing to do, just
#     rbind with fill so mismatched column sets don't silently drop rows)
#   - files with no ancestry column at all, in which case you supply
#     `ancestry_labels` aligned 1:1 with `files` and it's stamped on
# ==============================================================================
# ANCESTRY GROUP. This handles both flavors analysts send back:
#   - files that already contain an ancestry column (nothing to do, just
#     rbind with fill so mismatched column sets don't silently drop rows)
#   - files with no ancestry column at all, in which case you supply
#     `ancestry_labels` aligned 1:1 with `files` and it's stamped on
#
# skip_missing = FALSE (default) hard-fails if any file is missing, same as
# before. Set skip_missing = TRUE when you're templating a file list across
# many phenotypes/traits that don't all have every ancestry (e.g. building
# "${pheno}_eur.txt,${pheno}_afr.txt,${pheno}_eas.txt,..." in a loop and not
# every phenotype has an AFR or EAS file) - missing files (and their aligned
# label, if ancestry_labels was supplied) are dropped, with a warning naming
# exactly which ones were skipped, rather than the whole run aborting.
# ==============================================================================
combine_ancestry_files <- function(files, ancestry_col = "ancestry", ancestry_labels = NULL,
                                   skip_missing = FALSE) {
  files <- trimws(files)
  stopifnot(length(files) > 0)
  if (!is.null(ancestry_labels)) {
    ancestry_labels <- trimws(ancestry_labels)
    if (length(ancestry_labels) != length(files))
      stop("ancestry_labels must be the same length as files (one label per file), got ",
           length(ancestry_labels), " labels for ", length(files), " files.")
  }
  
  found <- file.exists(files)
  if (!all(found)) {
    missing_desc <- if (!is.null(ancestry_labels))
      paste0(files[!found], " (", ancestry_labels[!found], ")") else files[!found]
    if (!skip_missing)
      stop(sum(!found), " of ", length(files), " file(s) not found: ",
           paste(missing_desc, collapse = ", "),
           ". Pass skip_missing = TRUE (--skip-missing-ancestries on the CLI) to skip ",
           "missing ancestries instead of failing the whole run.")
    warning(sum(!found), " of ", length(files), " ancestry file(s) not found and skipped: ",
            paste(missing_desc, collapse = ", "))
    if (!is.null(ancestry_labels)) ancestry_labels <- ancestry_labels[found]
    files <- files[found]
    if (length(files) == 0)
      stop("None of the supplied files exist; nothing to combine.")
  }
  
  parts <- lapply(seq_along(files), function(i) {
    dt <- fread(files[i])
    if (ancestry_col %in% names(dt)) {
      if (!is.null(ancestry_labels) && !all(dt[[ancestry_col]] == ancestry_labels[i])) {
        warning(files[i], " already has an `", ancestry_col, "` column; ignoring the ",
                "supplied label '", ancestry_labels[i], "' and keeping the file's own values.")
      }
    } else {
      if (is.null(ancestry_labels))
        stop(files[i], " has no `", ancestry_col, "` column; supply ancestry_labels.")
      dt[[ancestry_col]] <- ancestry_labels[i]
    }
    dt
  })
  rbindlist(parts, fill = TRUE, use.names = TRUE)
}

# ==============================================================================
# recode_sex() - normalize a site-specific sex code into "Male"/"Female"
# character labels.
#
# Biobanks vary in how they code sex (0/1, 1/2, "M"/"F", ...), and that
# encoding can even differ across ancestry-specific files from the SAME
# biobank. Rather than have every downstream group_by() silently trust
# whatever numeric convention showed up, this maps explicit
# male_value/female_value onto "Male"/"Female" and turns anything else into
# NA (with a warning) instead of guessing.
#
# No-op (with no warning) if the column already holds "Male"/"Female"
# strings, so it's safe to call unconditionally even on pre-cleaned data.
# ==============================================================================
recode_sex <- function(df, sex_col = "sex", male_value = 1, female_value = 0) {
  raw <- df[[sex_col]]
  if (is.character(raw) && all(raw %in% c("Male", "Female", NA_character_))) return(df)
  
  recoded <- rep(NA_character_, length(raw))
  recoded[as.character(raw) == as.character(male_value)]   <- "Male"
  recoded[as.character(raw) == as.character(female_value)] <- "Female"
  
  n_unmapped <- sum(is.na(recoded) & !is.na(raw))
  if (n_unmapped > 0)
    warning(n_unmapped, " `", sex_col, "` value(s) matched neither male_value (",
            male_value, ") nor female_value (", female_value, ") and were set to NA. ",
            "Distinct raw values seen: ", paste(unique(raw), collapse = ", "))
  
  df[[sex_col]] <- recoded
  df
}

# ==============================================================================
# .year_from_date_or_year() - parses a column that may hold either full dates
# ("1990-05-01") or bare 4-digit years ("1990" / 1990) into a numeric year.
#
# as.Date() has no reliable way to parse a bare year: as.Date("1990") errors,
# and as.Date(1990) "succeeds" by treating 1990 as a day-count offset from
# 1970-01-01 (landing around 1975) - wrong, with no warning. Bare-year values
# are detected first (by pattern, not by trying and catching) and taken at
# face value; only the remainder is run through as.Date().
# ==============================================================================
.year_from_date_or_year <- function(x) {
  x_chr <- trimws(as.character(x))
  is_bare_year <- grepl("^[0-9]{4}$", x_chr)
  
  year <- rep(NA_real_, length(x))
  year[is_bare_year] <- as.numeric(x_chr[is_bare_year])
  
  needs_date <- !is_bare_year & !is.na(x_chr) & nzchar(x_chr)
  if (any(needs_date)) {
    parsed <- suppressWarnings(as.Date(x_chr[needs_date]))
    year[needs_date] <- as.numeric(format(parsed, "%Y"))
  }
  
  n_unparsed <- sum(is.na(year) & !is.na(x_chr) & nzchar(x_chr))
  if (n_unparsed > 0)
    warning(n_unparsed, " value(s) could not be parsed as a date (expects ",
            "YYYY-MM-DD) or a bare 4-digit year, and were set to NA. Example ",
            "unparsed value: '", x_chr[which(is.na(year) & !is.na(x_chr) & nzchar(x_chr))[1]], "'")
  year
}

# ==============================================================================
# .convert_to_years() - convert an age/time column that some biobanks record
# in days or months (instead of years) into years, at the point it first
# enters the script. Everything downstream (birthyear + age, Surv() times,
# month_to_progression = year_to_progression * 12, histogram bin widths,
# KM curve x-axes, ...) assumes years, so converting once here - rather than
# threading a units flag through every consumer - is what keeps that
# assumption valid without touching the rest of the script.
#
# `units` matches the CLI's --units flag: "years" (default, no-op),
# "days" (divide by 365), or "months" (divide by 12). Both divisors are
# fixed constants (not 365.25, not calendar-month-aware) - fine for the
# summary-statistics/histogram/KM use case here, but worth knowing if exact
# day/month arithmetic ever matters upstream.
# ==============================================================================
.convert_to_years <- function(x, units = c("years", "days", "months")) {
  units <- match.arg(units)
  divisor <- switch(units, years = 1, days = 365, months = 12)
  x / divisor
}

# ==============================================================================
# .join_age_from_file() - shared merge logic used by both join_enrollment()
# and join_followup(): read `file`, pull an age column or (birthyear-derived)
# date column, and left-join it onto `cohort` under `out_col`. Not exported;
# called by the two functions below so the date-vs-age handling only lives
# in one place.
# ==============================================================================
.join_age_from_file <- function(cohort, file, id_col, file_id_col,
                                age_col = NULL, date_col = NULL,
                                birthyear_col = "birthyear", out_col) {
  cohort <- as.data.table(cohort)
  src <- fread(file)
  if (!is.null(age_col)) {
    src <- src[, .(.JOIN_ID = get(file_id_col), .JOIN_AGE = get(age_col))]
  } else if (!is.null(date_col)) {
    join_year <- .year_from_date_or_year(src[[date_col]])
    src <- src[, .(.JOIN_ID = get(file_id_col), .JOIN_YEAR = join_year)]
  } else {
    stop("Supply exactly one of an age column or a date column from ", file, ".")
  }
  setnames(src, ".JOIN_ID", id_col)
  cohort <- merge(cohort, src, by = id_col, all.x = TRUE)
  if (!is.null(date_col)) {
    cohort[, (out_col) := .JOIN_YEAR - get(birthyear_col)]
    cohort[, .JOIN_YEAR := NULL]
  } else {
    cohort[, (out_col) := .JOIN_AGE]
    cohort[, .JOIN_AGE := NULL]
  }
  cohort
}

# ==============================================================================
# join_enrollment() - attach enrollment age from a separate file, since
# biobanks are not expected to include enrollment time in the phenotype file
# itself. Accepts either an enrollment AGE column directly, or an enrollment
# DATE (converted to age using birthyear).
# ==============================================================================
join_enrollment <- function(cohort,
                            enrollment_file,
                            id_col = "IID",
                            enrollment_id_col = "IID",
                            enrollment_age_col = NULL,      # age, if already an age
                            enrollment_date_col = NULL, # date, converted via birthyear
                            birthyear_col = "birthyear",
                            out_col = "enrollment_age") {
  cohort <- .join_age_from_file(cohort, enrollment_file, id_col, enrollment_id_col,
                                age_col = enrollment_age_col, date_col = enrollment_date_col,
                                birthyear_col = birthyear_col, out_col = out_col)
  n_missing <- sum(is.na(cohort[[out_col]]))
  if (n_missing > 0)
    message(n_missing, " of ", nrow(cohort), " people are missing `", out_col,
            "` after the join - check `", id_col, "` overlap between the cohort ",
            "file and ", enrollment_file)
  cohort[]
}

# ==============================================================================
# join_followup() - attach last-followup age for ALL individuals from a
# separate file. This may be the SAME file as enrollment (e.g. one
# demographics extract with both enrollment and last-contact/death dates) or
# a DIFFERENT one (e.g. a dedicated vital-status file) - pass whichever
# applies via followup_file; enrollment_file and followup_file are read and
# joined independently either way, so pointing both at the same path is
# fine.
#
# Handles the common GBMI phenotype convention where the pheno file's own
# age column is only a true "last followup" for controls (a case's recorded
# age is age-at-diagnosis instead). If `existing_col` already has a value
# for a person (e.g. carried over from the ancestry-specific phenotype
# file), the joined file's value is used only where it's non-missing;
# existing values are kept as a fallback rather than being wiped to NA. That
# means:
#   - a full-cohort followup file (covers everyone) simply replaces the
#     existing column outright;
#   - a cases-only or partial vital-status file fills in just those rows,
#     leaving the age-at-diagnosis proxy in place for whoever it doesn't
#     cover - with a warning so that's not silently invisible.
# ==============================================================================
join_followup <- function(cohort,
                          followup_file,
                          id_col = "IID",
                          followup_id_col = "IID",
                          followup_col = NULL,      # age, if already an age
                          followup_date_col = NULL, # date, converted via birthyear
                          birthyear_col = "birthyear",
                          existing_col = NULL,
                          out_col = "followup_age") {
  cohort <- as.data.table(cohort)
  had_existing <- !is.null(existing_col) && existing_col %in% names(cohort)
  cohort <- .join_age_from_file(cohort, followup_file, id_col, followup_id_col,
                                age_col = followup_col, date_col = followup_date_col,
                                birthyear_col = birthyear_col, out_col = ".JOINED_FOLLOWUP")
  
  if (had_existing) {
    n_filled_only <- sum(is.na(cohort[[".JOINED_FOLLOWUP"]]) & !is.na(cohort[[existing_col]]))
    if (n_filled_only > 0)
      warning(n_filled_only, " of ", nrow(cohort), " people have no match in ",
              followup_file, "; keeping the existing `", existing_col,
              "` value for those rows (this is a proxy - e.g. age-at-diagnosis - ",
              "unless you know otherwise).")
    cohort[, (out_col) := fifelse(!is.na(.JOINED_FOLLOWUP), .JOINED_FOLLOWUP, get(existing_col))]
  } else {
    n_missing <- sum(is.na(cohort[[".JOINED_FOLLOWUP"]]))
    if (n_missing > 0)
      message(n_missing, " of ", nrow(cohort), " people are missing `", out_col,
              "` after the join - check `", id_col, "` overlap between the cohort ",
              "file and ", followup_file)
    cohort[, (out_col) := .JOINED_FOLLOWUP]
  }
  cohort[, .JOINED_FOLLOWUP := NULL]
  cohort[]
}

# clean up KM curve in case of privacy issues
group_km <- function(onset_km, n_min = 5, 
                     group_cols = c("ancestry", "sex"), 
                     min_risk = n_min) {
  
  bin_one_group <- function(km_df) {
    km_df <- km_df %>% arrange(time)
    
    bin_id <- integer(nrow(km_df))
    current_bin <- 1L
    current_n <- 0L
    current_event <- 0L
    current_censor <- 0L
    
    for (i in seq_len(nrow(km_df))) {
      
      bin_id[i] <- current_bin
      
      current_n <- current_n + km_df$n.event[i] + km_df$n.censor[i]
      current_event <- current_event + km_df$n.event[i]
      current_censor <- current_censor + km_df$n.censor[i]
      
      if (current_n >= n_min & current_event >= n_min & current_censor >= n_min) {
        current_bin <- current_bin + 1L
        current_n <- 0L
        current_event <- 0L
        current_censor <- 0L
      }
    }
    
    bin_check <- tibble(
      bin = bin_id,
      n_event = km_df$n.event,
      n_censor = km_df$n.censor
    ) %>%
      group_by(bin) %>%
      summarise(
        n_bin = sum(n_event + n_censor, na.rm = TRUE),
        n_event = sum(n_event, na.rm = TRUE),
        n_censor = sum(n_censor, na.rm = TRUE),
        .groups = "drop"
      )
    
    if (nrow(bin_check) > 1) {
      last_bin <- max(bin_check$bin)
      last_row <- bin_check %>% filter(bin == last_bin)
      
      last_bin_fails <- with(
        last_row,
        n_bin < n_min | n_event < n_min | n_censor < n_min
      )
      
      if (isTRUE(last_bin_fails)) {
        bin_id[bin_id == last_bin] <- last_bin - 1L
      }
    }
    
    km_df %>%
      mutate(bin = bin_id) %>%
      group_by(bin) %>%
      summarise(
        time_start = first(time),
        time_end   = last(time),
        time       = last(time),
        
        n.risk   = first(n.risk),
        n.event  = sum(n.event),
        n.censor = sum(n.censor),
        
        estimate  = last(estimate),
        std.error = last(std.error),
        conf.low  = last(conf.low),
        conf.high = last(conf.high),
        
        n_bin = sum(n.event + n.censor),
        
        across(
          any_of(c("n_total", "n_events_total")),
          first
        ),
        
        .groups = "drop"
      ) %>%
      filter(
        n_bin >= n_min,
        n.risk >= min_risk
      )
  }
  
  onset_km %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(~ bin_one_group(.x)) %>%
    ungroup()
}

# ==============================================================================
# Shared "combine sex/ancestry into ALL rows" step, factored out since it was
# duplicated three times in the original script (and the original bug - all
# three copies checked `cohort_clean$ancestry` instead of the local table -
# is fixed here by taking the table as an explicit argument).
# ==============================================================================
add_all_rows <- function(df, ancestry_col = "ancestry") {
  if (length(unique(df[[ancestry_col]])) == 1) {
    rbind(df, df %>% mutate(sex = "ALL"))
  } else {
    rbind(df,
          df %>% mutate(sex = "ALL"),
          df %>% mutate(ancestry = "ALL"),
          df %>% mutate(sex = "ALL", ancestry = "ALL"))
  }
}

# Apply consortium disclosure suppression to a long "group -> n" histogram
# table, replacing the original's paste0('<', count_min) string hack (which
# turned `n` into a character column and would silently break any downstream
# numeric use of it). Uses the same NA + minimum_cell_count_pass convention
# as the rest of the pipeline.
suppress_histogram <- function(df, min_cell_count) {
  if (min_cell_count <= 0) return(df)
  as.data.table(df) %>%
    suppress_small_cells(count_cols = "n", threshold = min_cell_count) %>%
    as_tibble()
}

#####################
## OVERALL SUMMARY ##
#####################
# Builds the combined cohort dataset (ancestry, sex, enrollment age, follow-up
# age) and writes the summary/histogram tables. Pass either a single
# already-combined `cohort` data.frame, OR `pheno_files` (+ optional
# `ancestry_labels`) to have combine_ancestry_files() build it for you.
#
# `enrollment_file` / `followup_file` may be the same path or different paths
# - each is joined independently, so a single demographics extract that has
# both enrollment and last-contact columns works by pointing both arguments
# at it with the right column names.
overall_summary <- function(cohort = NULL,
                            prefix,
                            pheno_files = NULL,
                            ancestry_labels = NULL,
                            enrollment_file = NULL,
                            enrollment_id_col = "IID",
                            # enrollment_col = NULL,
                            enrollment_date_col = NULL,
                            followup_file = NULL,
                            followup_id_col = "IID",
                            followup_col_in_file = NULL,
                            followup_date_col_in_file = NULL,
                            followup_date_col = NULL,
                            id_col = 'IID',
                            sex_col = 'sex',
                            ancestry_col = 'ancestry',
                            birthyear_col = 'birthyear',
                            units = 'years',
                            male_value = 1,
                            female_value = 0,
                            enrollment_age_col = NULL,
                            followup_col = 'age_death_or_lastvisit',
                            skip_missing_ancestries = FALSE,
                            min_cell_count = 0) {
  
  units <- match.arg(units, c("years", "days", "months"))
  
  if (is.null(cohort)) {
    stopifnot(!is.null(pheno_files))
    cohort <- combine_ancestry_files(pheno_files, ancestry_col = ancestry_col,
                                     ancestry_labels = ancestry_labels,
                                     skip_missing = skip_missing_ancestries)
  }
  if (!is.null(enrollment_file)) {
    cohort <- join_enrollment(cohort, enrollment_file,
                              id_col = id_col, enrollment_id_col = enrollment_id_col,
                              enrollment_age_col = enrollment_age_col,
                              enrollment_date_col = enrollment_date_col,
                              birthyear_col = birthyear_col,
                              out_col = "enrollment_age")
    enrollment_age_col <- "enrollment_age"
  } else if (!is.null(enrollment_date_col)) {
    # No separate file to join - the date (or bare year) column is already in
    # the combined pheno/cohort table, so convert it to age in place using
    # birthyear rather than routing through join_enrollment()/a file read.
    cohort <- as.data.table(cohort)
    cohort[, enrollment_age := .year_from_date_or_year(get(enrollment_date_col)) - get(birthyear_col)]
    enrollment_age_col <- "enrollment_age"
  }
  if (!is.null(followup_file)) {
    cohort <- join_followup(cohort, followup_file,
                            id_col = id_col, followup_id_col = followup_id_col,
                            followup_col = followup_col_in_file,
                            followup_date_col = followup_date_col_in_file,
                            birthyear_col = birthyear_col,
                            existing_col = followup_col,
                            out_col = followup_col)
  } else if (!is.null(followup_date_col)) {
    # Same idea for follow-up: a date (or bare year) column already in the
    # combined table, converted to age in place - no file to join.
    cohort <- as.data.table(cohort)
    cohort[, (followup_col) := .year_from_date_or_year(get(followup_date_col)) - get(birthyear_col)]
  }
  
  # rename columns
  cohort_clean <- cohort %>%
    select(all_of(c(id_col, sex_col, ancestry_col,
                    birthyear_col, enrollment_age_col, followup_col))) %>%
    rename(ID = all_of(id_col),
           sex = all_of(sex_col),
           ancestry = all_of(ancestry_col),
           birthyear = all_of(birthyear_col),
           enrollment_age = all_of(enrollment_age_col),
           followup_age = all_of(followup_col)) %>%
    mutate(followup_time = followup_age - enrollment_age,
           enrollment_year = birthyear + enrollment_age)
  
  cohort_clean <- recode_sex(cohort_clean, "sex", male_value, female_value)
  
  # sample sizes
  # cohort_n <- cohort_clean %>%
  #   group_by(ancestry, sex) %>%
  #   summarize(n = n(), .groups = 'drop') %>%
  #   suppress_histogram(min_cell_count)
  
  cohort_combined <- add_all_rows(cohort_clean, "ancestry")
  
  cohort_n <- cohort_combined %>%
    group_by(ancestry, sex) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  cohort_summary <- cohort_combined %>%
    group_by(ancestry, sex) %>%
    summarize_at(.vars = vars(birthyear, enrollment_age, enrollment_year, followup_age, followup_time),
                 .funs = list(mean = 'mean', median = 'median', sd = 'sd')) %>%
    ungroup()
  
  # histograms
  birthyear_histogram <- cohort_combined %>%
    mutate(birthyear = round(birthyear, 0)) %>%
    group_by(ancestry, sex, birthyear) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  followup_histogram <- cohort_combined %>%
    mutate(followup = round(followup_time)) %>%
    group_by(ancestry, sex, followup) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  enrollmentYear_histogram <- cohort_combined %>%
    mutate(enrollment_year = round(enrollment_year)) %>%
    group_by(ancestry, sex, enrollment_year) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  enrollmentAge_histogram <- cohort_combined %>%
    mutate(enrollment_age = round(enrollment_age)) %>%
    group_by(ancestry, sex, enrollment_age) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  fwrite(cohort_n, paste0(prefix, '_sample_sizes.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(cohort_summary, paste0(prefix, '_summaries.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(birthyear_histogram, paste0(prefix, '_birthyear_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(followup_histogram, paste0(prefix, '_followup_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(enrollmentYear_histogram, paste0(prefix, '_enrollmentYear_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(enrollmentAge_histogram, paste0(prefix, '_enrollmentAge_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  invisible(cohort_clean)
}

###########
## ONSET ##
###########
# Pass either `phenofile` (already combined, with an ancestry column) or
# `pheno_files` (+ optional `ancestry_labels`) to have combine_ancestry_files()
# stack the ancestry-specific files for you.
onset_summary <- function(phenofile = NULL,
                          prefix,
                          pheno_files = NULL,
                          ancestry_labels = NULL,
                          id_col = 'IID',
                          event_col = 'event',
                          eventage_col = 'age', # age of diagnosis for cases, age of last followup for controls
                          birthyear_col = 'birthyear',
                          sex_col = 'sex',
                          ancestry_col = 'ancestry',
                          male_value = 1,
                          female_value = 0,
                          units = 'years', # units of eventage_col in the input file(s): years (default), days, or months
                          km_min = 0, km_round = F,
                          min_cell_count = 0,
                          skip_missing_ancestries = FALSE) {
  
  units <- match.arg(units, c("years", "days", "months"))
  
  if (is.null(phenofile)) {
    stopifnot(!is.null(pheno_files))
    phenofile <- combine_ancestry_files(pheno_files, ancestry_col = ancestry_col,
                                        ancestry_labels = ancestry_labels,
                                        skip_missing = skip_missing_ancestries)
  }
  
  phenofile_clean <- phenofile %>%
    select(all_of(c(id_col, sex_col, ancestry_col, event_col, birthyear_col, eventage_col))) %>%
    rename(ID = all_of(id_col),
           sex = all_of(sex_col),
           event = all_of(event_col),
           ancestry = all_of(ancestry_col),
           birthyear = all_of(birthyear_col),
           age = all_of(eventage_col)) %>%
    mutate(age = .convert_to_years(age, units),
           eventyear = birthyear + age)
  
  phenofile_clean <- recode_sex(phenofile_clean, "sex", male_value, female_value)
  
  # onset_n <- phenofile_clean %>%
  #   group_by(ancestry, sex, event) %>%
  #   summarize(n = n(), .groups = 'drop') %>%
  #   suppress_histogram(min_cell_count)
  
  onset_combined <- add_all_rows(phenofile_clean, "ancestry")
  
  onset_n <- onset_combined %>%
    group_by(ancestry, sex, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  onset_summary_tbl <- onset_combined %>%
    group_by(ancestry, sex, event) %>%
    summarize_at(.vars = vars(age, eventyear),
                 .funs = list(mean = 'mean', median = 'median', sd = 'sd')) %>%
    ungroup()
  
  onset_age_histogram <- onset_combined %>%
    mutate(age = round(age)) %>%
    group_by(ancestry, sex, age, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  onset_year_histogram <- onset_combined %>%
    mutate(eventyear = round(eventyear)) %>%
    group_by(ancestry, sex, eventyear, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  if (km_round) {
    onset_km <- onset_combined %>%
      mutate(age = round(age)) %>%
      group_by(ancestry, sex) %>%
      nest() %>%
      mutate(
        fit = map(data, ~ survfit(Surv(age, event) ~ 1, data = .x)),
        km = map(fit, broom::tidy),
        n_total = map_int(data, nrow),
        n_events_total = map_int(data, ~ sum(.x$event == 1, na.rm = TRUE))
      ) %>%
      select(-data, -fit) %>%
      unnest(km) %>%
      ungroup()
  } else {
    onset_km <- onset_combined %>%
      group_by(ancestry, sex) %>%
      nest() %>%
      mutate(
        fit = map(data, ~ survfit(Surv(age, event) ~ 1, data = .x)),
        km = map(fit, broom::tidy),
        n_total = map_int(data, nrow),
        n_events_total = map_int(data, ~ sum(.x$event == 1, na.rm = TRUE))
      ) %>%
      select(-data, -fit) %>%
      unnest(km) %>%
      ungroup()
  }
  
  if (km_min > 0) onset_km <- group_km(onset_km, n_min = km_min)
  
  fwrite(onset_n, paste0(prefix, '_sample_sizes.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(onset_summary_tbl, paste0(prefix, '_summaries.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(onset_age_histogram, paste0(prefix, '_age_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(onset_year_histogram, paste0(prefix, '_year_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(onset_km, paste0(prefix, '_KM.txt.gz'), sep = '\t', compress = 'gzip')
  invisible(phenofile_clean)
}

#################
## PROGRESSION ##
#################
prog_summary <- function(phenofile = NULL,
                         prefix,
                         pheno_files = NULL,
                         ancestry_labels = NULL,
                         id_col = 'IID',
                         event_col = 'secondEvent',
                         eventage_col = 'ageExit', # age of T2 for cases, age of last followup for controls
                         time_to_progression_col = 'secondTime', # fallback: onsetAge + this, used only if eventage_col isn't in the file
                         onsetage_col = 'diagAge', # age of T1 diagnosis
                         birthyear_col = 'birthyear',
                         sex_col = 'sex',
                         ancestry_col = 'ancestry',
                         male_value = 1,
                         female_value = 0,
                         units = 'years,years', # units of eventage_col/onsetage_col/time_to_progression_col in the input file(s): years (default), days, or months
                         km_min = 0, km_round = F,
                         min_cell_count = 0,
                         skip_missing_ancestries = FALSE) {
  
  unit_split = strsplit(units, split=',')[[1]]
  units <- unname(sapply(unit_split, FUN=function(x) match.arg(x, c("years", "days", "months"))))
  
  if (is.null(phenofile)) {
    stopifnot(!is.null(pheno_files))
    phenofile <- combine_ancestry_files(pheno_files, ancestry_col = ancestry_col,
                                        ancestry_labels = ancestry_labels,
                                        skip_missing = skip_missing_ancestries)
  }
  
  # Both the baseline (onsetage_col) and event (eventage_col, or the
  # time_to_progression_col fallback below) time columns are converted to
  # years here, at the point they first enter the script, using the same
  # `units` the rest of the file's columns are in - onsetage_col is always
  # present so it's converted unconditionally, before eventage_col/
  # time_to_progression_col is even resolved, so the fallback derivation
  # below (onsetAge + time-to-progression) is never a mismatched-units sum.
  phenofile <- as.data.table(phenofile)
  phenofile[, (onsetage_col) := .convert_to_years(get(onsetage_col), units[1])]
  
  # eventage_col (an absolute age at the second event) is preferred when the
  # phenotype file actually has it. If it doesn't - some sites only record
  # time-to-progression (e.g. `secondTime`) rather than an age - fall back
  # to deriving it as onsetAge + time_to_progression_col instead.
  if (!is.null(eventage_col) && eventage_col %in% names(phenofile)) {
    phenofile[, (eventage_col) := .convert_to_years(get(eventage_col), units[2])]
  } else if (!is.null(time_to_progression_col) && time_to_progression_col %in% names(phenofile)) {
    if (!is.null(eventage_col))
      message("`", eventage_col, "` not found in the phenotype file(s); deriving age at ",
              "second event as `", onsetage_col, "` + `", time_to_progression_col, "` instead.")
    phenofile[, .DERIVED_EVENTAGE := get(onsetage_col) + .convert_to_years(get(time_to_progression_col), units[2])]
    eventage_col <- ".DERIVED_EVENTAGE"
  } else {
    stop("Neither eventage_col ('", eventage_col, "') nor time_to_progression_col ('",
         time_to_progression_col, "') is a column in the phenotype file(s). Supply ",
         "one via --eventage-col (an absolute age) or --time-to-progression-col ",
         "(elapsed time since onset, in the units given by --units, added to ",
         "--onsetage-col).")
  }
  
  phenofile_clean <- phenofile %>%
    select(all_of(c(id_col, sex_col, ancestry_col, event_col,
                    eventage_col, onsetage_col, birthyear_col))) %>%
    rename(ID = all_of(id_col),
           sex = all_of(sex_col),
           ancestry = all_of(ancestry_col),
           event = all_of(event_col),
           eventAge = all_of(eventage_col),
           onsetAge = all_of(onsetage_col),
           birthyear = all_of(birthyear_col)) %>%
    mutate(eventyear = birthyear + eventAge,
           year_to_progression = eventAge - onsetAge,
           month_to_progression = year_to_progression * 12)
  
  phenofile_clean <- recode_sex(phenofile_clean, "sex", male_value, female_value)
  
  # prog_n <- phenofile_clean %>%
  #   group_by(ancestry, sex, event) %>%
  #   summarize(n = n(), .groups = 'drop') %>%
  #   suppress_histogram(min_cell_count)
  
  prog_combined <- add_all_rows(phenofile_clean, "ancestry")
  
  prog_n <- prog_combined %>%
    group_by(ancestry, sex, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  prog_summary_tbl <- prog_combined %>%
    group_by(ancestry, sex, event) %>%
    summarize_at(.vars = vars(onsetAge, eventAge, eventyear),
                 .funs = list(mean = 'mean', median = 'median', sd = 'sd')) %>%
    ungroup()
  
  prog_baseline_histogram <- prog_combined %>%
    mutate(baseline_age = round(onsetAge)) %>%
    group_by(ancestry, sex, baseline_age, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  prog_age_histogram <- prog_combined %>%
    mutate(event_age = round(eventAge)) %>%
    group_by(ancestry, sex, event_age, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  prog_year_histogram <- prog_combined %>%
    mutate(year = round(year_to_progression)) %>%
    group_by(ancestry, sex, year, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  prog_month_histogram <- prog_combined %>%
    mutate(month = round(month_to_progression)) %>%
    group_by(ancestry, sex, month, event) %>%
    summarize(n = n(), .groups = 'drop') %>%
    suppress_histogram(min_cell_count)
  
  fit_km <- function(time_col, km_round) {
    if (km_round) {
      prog_combined %>%
        mutate_at(.vars=vars(any_of(time_col)), .funs='round') %>%
        group_by(ancestry, sex) %>%
        nest() %>%
        mutate(
          fit = map(data, ~ survfit(as.formula(paste0("Surv(", time_col, ", event) ~ 1")), data = .x)),
          km = map(fit, broom::tidy),
          n_total = map_int(data, nrow),
          n_events_total = map_int(data, ~ sum(.x$event == 1, na.rm = TRUE))
        ) %>%
        select(-data, -fit) %>%
        unnest(km) %>%
        ungroup()
    } else {
      prog_combined %>%
        group_by(ancestry, sex) %>%
        nest() %>%
        mutate(
          fit = map(data, ~ survfit(as.formula(paste0("Surv(", time_col, ", event) ~ 1")), data = .x)),
          km = map(fit, broom::tidy),
          n_total = map_int(data, nrow),
          n_events_total = map_int(data, ~ sum(.x$event == 1, na.rm = TRUE))
        ) %>%
        select(-data, -fit) %>%
        unnest(km) %>%
        ungroup()
    }
    
  }
  prog_age_km <- fit_km("eventAge", km_round)
  prog_year_km <- fit_km("year_to_progression", km_round)
  prog_month_km <- fit_km("month_to_progression", km_round)
  
  if (km_min > 0) {
    prog_age_km <- group_km(prog_age_km, n_min = km_min)
    prog_year_km <- group_km(prog_year_km, n_min = km_min)
    prog_month_km <- group_km(prog_month_km, n_min = km_min)
  }
  
  fwrite(prog_n, paste0(prefix, '_sample_sizes.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_summary_tbl, paste0(prefix, '_summaries.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_baseline_histogram, paste0(prefix, '_baselineAge_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_age_histogram, paste0(prefix, '_eventAge_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_year_histogram, paste0(prefix, '_year_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_month_histogram, paste0(prefix, '_month_histogram.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_age_km, paste0(prefix, '_age_KM.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_year_km, paste0(prefix, '_year_KM.txt.gz'), sep = '\t', compress = 'gzip')
  fwrite(prog_month_km, paste0(prefix, '_month_KM.txt.gz'), sep = '\t', compress = 'gzip')
  invisible(phenofile_clean)
}

# ==============================================================================
# COMMAND-LINE ENTRY POINT
# ==============================================================================
is_main <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  any(grepl("^--file=", args)) && !interactive()
}

if (is_main()) {
  option_list <- list(
    make_option("--mode", type = "character",
                help = "One of: cohort, onset, progression [REQUIRED]"),
    make_option("--pheno-files", type = "character",
                help = "Comma-separated ancestry-specific input files [REQUIRED]"),
    make_option("--ancestry-labels", type = "character", default = NULL,
                help = "Comma-separated labels aligned with --pheno-files; omit if files already have an ancestry column"),
    make_option("--skip-missing-ancestries", action = "store_true", default = FALSE,
                help = "Skip --pheno-files entries that don't exist instead of failing the run - useful when templating a file list (e.g. ${pheno}_eur.txt,${pheno}_afr.txt,...) across phenotypes that don't all have every ancestry"),
    make_option("--km-round", action = "store_true", default = FALSE,
                help = "round KM time points to nearest integer"),
    make_option("--prefix", type = "character", help = "Output file prefix [REQUIRED]"),
    make_option("--min-cell-count", type = "integer", default = 0,
                help = "Disclosure: suppress counts below this (default: 0)"),
    make_option("--km-min", type = "integer", default = 0,
                help = "Minimum bin size for KM curve grouping, onset/progression only (default: 0, off)"),
    make_option("--units", type = "character", default = NULL,
                help = "onset/progression only: units the event-time column(s) are recorded in the input file(s) - one of years (default), days, or months. Applies to --eventage-col for onset; to --eventage-col, --onsetage-col, and --time-to-progression-col for progression. Converted to years internally (days / 365, months / 12) before any downstream calculation."),
    make_option("--male-value", type = "character", default = "1",
                help = "Raw sex code meaning Male in the input file(s) (default: 1)"),
    make_option("--female-value", type = "character", default = "0",
                help = "Raw sex code meaning Female in the input file(s) (default: 0)"),
    # cohort mode only
    make_option("--enrollment-file", type = "character", default = NULL,
                help = "cohort mode: file with enrollment age or date"),
    make_option("--enrollment-id-col", type = "character", default = "IID"),
    make_option("--enrollment-col", type = "character", default = NULL,
                help = "cohort mode: enrollment AGE column name in --enrollment-file"),
    make_option("--enrollment-date-col", type = "character", default = NULL,
                help = "cohort mode: enrollment DATE column name in --enrollment-file (converted via birthyear)"),
    make_option("--followup-file", type = "character", default = NULL,
                help = "cohort mode: file with last-followup age or date for some/all individuals; may be the same file as --enrollment-file"),
    make_option("--followup-id-col", type = "character", default = "IID"),
    make_option("--followup-col-in-file", type = "character", default = NULL,
                help = "cohort mode: last-followup AGE column name in --followup-file"),
    make_option("--followup-date-col-in-file", type = "character", default = NULL,
                help = "cohort mode: last-followup DATE column name in --followup-file (converted via birthyear)"),
    make_option("--followup-date-col", type = "character", default = NULL,
                help = "cohort mode: last-followup DATE (or bare-year) column name already in --pheno-files, used only when --followup-file is NOT supplied"),
    # shared column-name overrides (defaults match the GBMI phenotype convention)
    make_option("--id-col", type = "character", default = "IID"),
    make_option("--sex-col", type = "character", default = "sex"),
    make_option("--ancestry-col", type = "character", default = "ancestry"),
    make_option("--birthyear-col", type = "character", default = "birthyear"),
    make_option("--followup-col", type = "character", default = "age_death_or_lastvisit"),
    make_option("--event-col", type = "character", default = "event"),
    make_option("--eventage-col", type = "character", default = "age"),
    make_option("--onsetage-col", type = "character", default = "diagAge"),
    make_option("--time-to-progression-col", type = "character", default = "secondTime",
                help = "progression mode: fallback used only if --eventage-col isn't found in the phenotype file(s); eventAge is derived as onsetAge + this (years)")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$mode) || is.null(opt[["pheno-files"]]) || is.null(opt$prefix))
    stop("--mode, --pheno-files, and --prefix are required.")
  
  pheno_files <- strsplit(opt[["pheno-files"]], ",")[[1]]
  ancestry_labels <- if (!is.null(opt[["ancestry-labels"]])) strsplit(opt[["ancestry-labels"]], ",")[[1]] else NULL
  
  if (is.null(opt$units)) {
    opt$units = ifelse(opt$mode=='progression', 'years,years', 'years')
  }
  if (!opt$units %in% c("years", "days", "months", 
                        apply(expand.grid(x1 = c("years", "days", "months"), x2 = c("years", "days", "months")), 
                              1, paste, collapse=',')))
    stop("--units must be one of: years, days, months (got '", opt$units, "')")
  
  # create directory as needed
  if (!dir.exists(dirname(opt$prefix))) {
    cat('creating output directory ', dirname(opt$prefix), '\n')
    dir.create(dirname(opt$prefix), recursive = T)
  }
  
  if (opt$mode == "cohort") {
    overall_summary(prefix = opt$prefix, pheno_files = pheno_files, ancestry_labels = ancestry_labels,
                    enrollment_file = opt[["enrollment-file"]], enrollment_id_col = opt[["enrollment-id-col"]],
                    enrollment_age_col = opt[["enrollment-col"]], enrollment_date_col = opt[["enrollment-date-col"]],
                    followup_file = opt[["followup-file"]], followup_id_col = opt[["followup-id-col"]],
                    followup_col_in_file = opt[["followup-col-in-file"]],
                    followup_date_col_in_file = opt[["followup-date-col-in-file"]],
                    followup_date_col = opt[["followup-date-col"]],
                    id_col = opt[["id-col"]], sex_col = opt[["sex-col"]], ancestry_col = opt[["ancestry-col"]],
                    birthyear_col = opt[["birthyear-col"]],
                    units = opt[["units"]],
                    male_value = opt[["male-value"]], female_value = opt[["female-value"]],
                    followup_col = opt[["followup-col"]], min_cell_count = opt[["min-cell-count"]],
                    skip_missing_ancestries = opt[["skip-missing-ancestries"]])
  } else if (opt$mode == "onset") {
    onset_summary(prefix = opt$prefix, pheno_files = pheno_files, ancestry_labels = ancestry_labels,
                  id_col = opt[["id-col"]], event_col = opt[["event-col"]], eventage_col = opt[["eventage-col"]],
                  birthyear_col = opt[["birthyear-col"]], sex_col = opt[["sex-col"]], ancestry_col = opt[["ancestry-col"]],
                  male_value = opt[["male-value"]], female_value = opt[["female-value"]],
                  units = opt[["units"]],
                  km_min = opt[["km-min"]], km_round = opt[["km-round"]], min_cell_count = opt[["min-cell-count"]],
                  skip_missing_ancestries = opt[["skip-missing-ancestries"]])
  } else if (opt$mode == "progression") {
    prog_summary(prefix = opt$prefix, pheno_files = pheno_files, ancestry_labels = ancestry_labels,
                 id_col = opt[["id-col"]], event_col = opt[["event-col"]], eventage_col = opt[["eventage-col"]],
                 time_to_progression_col = opt[["time-to-progression-col"]],
                 onsetage_col = opt[["onsetage-col"]], birthyear_col = opt[["birthyear-col"]], sex_col = opt[["sex-col"]],
                 ancestry_col = opt[["ancestry-col"]],
                 male_value = opt[["male-value"]], female_value = opt[["female-value"]],
                 units = opt[["units"]],
                 km_min = opt[["km-min"]], km_round = opt[["km-round"]], min_cell_count = opt[["min-cell-count"]],
                 skip_missing_ancestries = opt[["skip-missing-ancestries"]])
  } else {
    stop("--mode must be one of: cohort, onset, progression (got '", opt$mode, "')")
  }
  cat("RUN COMPLETE: yes\n")
}
