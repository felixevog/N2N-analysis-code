rm(list = ls())   # clear environment so all objects are rebuilt from source files

# ===========================================================================
# N2N COMBINED DASS-21 EXTERNAL VALIDITY ANALYSIS
# Pools pilot + main study participants
#
# Three analyses (each run for unregulated and regulated PECAN):
#
#   Analysis 1 — pilot nar  +  main nar_ex  (explicit connections only)
#   Analysis 2 — pilot nar  +  main nar     (full: explicit + implicit)
#
# In both analyses the nar metrics for the pilot are identical (pilot has no
# explicit/implicit split). The distinction is which main-study nar variant
# is used for the 50 main-study participants.
#
# Common metrics correlated with DASS T1 and T2:
#   nar:   n_edges, strength, density, n_nodes, n_feedback_loops
#   PECAN: n_edges, strength, density, n_nodes, mean_severity,
#          n_feedback_loops (pecan_le4 unreg / pecan_reg_le4 reg)
#
# Output: output/combined/tables/
# ===========================================================================

library(tidyverse)

# ---------------------------------------------------------------------------
# 0. Load per-participant data from both studies
# ---------------------------------------------------------------------------

# Regulation rule (mirroring the study scripts):
# Keep at most floor(n_nodes × 1.0) edges per participant, strongest first;
# ties at the cut-off boundary are always kept in full.
regulate_pecan_edges <- function(edges_df, nodes_df, reg_value = 1.0) {
  n_nodes <- nrow(nodes_df)
  n_keep  <- floor(n_nodes * reg_value)
  if (nrow(edges_df) == 0 || n_keep <= 0) return(edges_df[integer(0), ])
  edges_sorted <- edges_df |> arrange(desc(width))
  if (nrow(edges_sorted) <= n_keep) return(edges_sorted)
  cutoff_width <- edges_sorted$width[n_keep]
  edges_sorted |> filter(width >= cutoff_width)
}

# -- 0a. Network metrics from long-format output files (pivoted to wide) --

pilot_long <- read.csv("output/pilotstudy/tables/4a_network_metrics_long.csv")
main_long  <- read.csv("output/mainstudy/tables/4a_ext_network_metrics_long.csv")

pilot_wide <- pilot_long |>
  pivot_wider(id_cols = vp, names_from = network,
              values_from = c(n_nodes, n_edges, density)) |>
  rename(
    nar_nodes         = n_nodes_nar,
    nar_n_edges       = n_edges_nar,
    nar_density       = density_nar,
    pecan_nodes       = n_nodes_PECAN_unreg,
    pecan_n_edges     = n_edges_PECAN_unreg,
    pecan_density     = density_PECAN_unreg,
    pecan_reg_n_edges = n_edges_PECAN_reg,
    pecan_reg_density = density_PECAN_reg
  ) |>
  mutate(nar_strength = nar_n_edges)  # nar edges are binary; strength = n_edges

main_wide <- main_long |>
  pivot_wider(id_cols = vp, names_from = network,
              values_from = c(n_nodes, n_edges, density)) |>
  rename(
    nar_nodes         = n_nodes_nar,
    nar_n_edges       = n_edges_nar,
    nar_density       = density_nar,
    nar_ex_nodes      = n_nodes_nar_ex,
    nar_ex_n_edges    = n_edges_nar_ex,
    nar_ex_density    = density_nar_ex,
    pecan_nodes       = n_nodes_PECAN_unreg,
    pecan_n_edges     = n_edges_PECAN_unreg,
    pecan_density     = density_PECAN_unreg,
    pecan_reg_n_edges = n_edges_PECAN_reg,
    pecan_reg_density = density_PECAN_reg
  ) |>
  mutate(nar_strength    = nar_n_edges,     # binary edges
         nar_ex_strength = nar_ex_n_edges)

pilot_vps <- pilot_wide$vp
main_vps  <- main_wide$vp

# -- 0b. PECAN strength (sum of edge widths) from raw edge files --
# nar edges are binary so nar_strength = nar_n_edges (handled above).
# PECAN edges carry a 1–5 width; regulated strength uses the filtered edge set.

stre_pilot_unreg <- map_dfr(pilot_vps, function(vp) {
  edges <- read.csv(paste0("data/bearbeitet/PECAN/", vp, "_PECAN_edges.csv"))
  tibble(vp = vp, pecan_strength = sum(edges$width, na.rm = TRUE))
})

stre_pilot_reg <- map_dfr(pilot_vps, function(vp) {
  edges <- read.csv(paste0("data/bearbeitet/PECAN/", vp, "_PECAN_edges.csv"))
  nodes <- read.csv(paste0("data/bearbeitet/PECAN/", vp, "_PECAN_nodes.csv"))
  tibble(vp = vp,
         pecan_strength = sum(regulate_pecan_edges(edges, nodes)$width, na.rm = TRUE))
})

stre_main_unreg <- map_dfr(main_vps, function(vp) {
  edges <- read.csv(paste0("data/bearbeitet/mainstudy/PECAN/", vp, "_PECAN_edges.csv"))
  tibble(vp = vp, pecan_strength = sum(edges$width, na.rm = TRUE))
})

stre_main_reg <- map_dfr(main_vps, function(vp) {
  edges <- read.csv(paste0("data/bearbeitet/mainstudy/PECAN/", vp, "_PECAN_edges.csv"))
  nodes <- read.csv(paste0("data/bearbeitet/mainstudy/PECAN/", vp, "_PECAN_nodes.csv"))
  tibble(vp = vp,
         pecan_strength = sum(regulate_pecan_edges(edges, nodes)$width, na.rm = TRUE))
})

# -- 0c. DASS scores --
# T1 files use VP01_01 as the VP column and DASS_total as the score column.

read_dass_t1 <- function(path) {
  read.csv(path) |> rename(vp = VP01_01, DASS_T1_total = DASS_total)
}

dass_pilot_t1 <- read_dass_t1("output/pilotstudy/tables/2a_DASS21_scores.csv")
dass_pilot_t2 <- read.csv("output/pilotstudy/tables/2a_DASS21_T2_scores.csv")
dass_main_t1  <- read_dass_t1("output/mainstudy/tables/2a_DASS21_scores.csv")
dass_main_t2  <- read.csv("output/mainstudy/tables/2a_DASS21_T2_scores.csv")

# -- 0d. Mean PECAN node severity per participant (from raw node files) --

sev_pilot <- map_dfr(pilot_vps, function(vp) {
  nodes <- read.csv(paste0("data/bearbeitet/PECAN/", vp, "_PECAN_nodes.csv"))
  tibble(vp = vp, pecan_mean_severity = mean(nodes$value, na.rm = TRUE))
})

sev_main <- map_dfr(main_vps, function(vp) {
  nodes <- read.csv(paste0("data/bearbeitet/mainstudy/PECAN/", vp, "_PECAN_nodes.csv"))
  tibble(vp = vp, pecan_mean_severity = mean(nodes$value, na.rm = TRUE))
})

# -- 0e. Feedback loops per participant --

loops_pilot <- read.csv("output/pilotstudy/tables/4d_feedback_loops.csv") |>
  dplyr::select(vp, pecan_le4, pecan_reg_le4, nar_le4)

loops_main <- read.csv("output/mainstudy/tables/4d_feedback_loops.csv") |>
  dplyr::select(vp, pecan_le4, pecan_reg_le4, nar_le4, nar_ex_le4)

# -- 0f. Assemble per-participant metric data frames --
# PECAN columns reflect unregulated or regulated edges depending on the variant.
#needed for partial correlations later on (weighted strength unnecessary)

nm_pilot_unreg <- pilot_wide |>
  transmute(vp, nar_n_edges, nar_strength, nar_density,
            pecan_n_edges, pecan_density) |>
  left_join(stre_pilot_unreg, by = "vp") |>
  left_join(dass_pilot_t1,    by = "vp") |>
  left_join(dass_pilot_t2,    by = "vp") |>
  mutate(study = "pilot")

nm_pilot_reg <- pilot_wide |>
  transmute(vp, nar_n_edges, nar_strength, nar_density,
            pecan_n_edges = pecan_reg_n_edges,
            pecan_density = pecan_reg_density) |>
  left_join(stre_pilot_reg, by = "vp") |>
  left_join(dass_pilot_t1,  by = "vp") |>
  left_join(dass_pilot_t2,  by = "vp") |>
  mutate(study = "pilot")

nm_main_unreg <- main_wide |>
  transmute(vp, nar_n_edges, nar_strength, nar_density,
            nar_ex_n_edges, nar_ex_strength, nar_ex_density,
            pecan_n_edges, pecan_density) |>
  left_join(stre_main_unreg, by = "vp") |>
  left_join(dass_main_t1,    by = "vp") |>
  left_join(dass_main_t2,    by = "vp") |>
  mutate(study = "main")

nm_main_reg <- main_wide |>
  transmute(vp, nar_n_edges, nar_strength, nar_density,
            nar_ex_n_edges, nar_ex_strength, nar_ex_density,
            pecan_n_edges = pecan_reg_n_edges,
            pecan_density = pecan_reg_density) |>
  left_join(stre_main_reg, by = "vp") |>
  left_join(dass_main_t1,  by = "vp") |>
  left_join(dass_main_t2,  by = "vp") |>
  mutate(study = "main")

# Node count tables (used by build_combined)
nc_pilot <- pilot_wide |> dplyr::select(vp, nar_nodes, pecan_nodes)
nc_main  <- main_wide  |> dplyr::select(vp, nar_nodes, nar_ex_nodes, pecan_nodes)

# ---------------------------------------------------------------------------
# 1. Build pooled datasets
# ---------------------------------------------------------------------------

# DASS + PECAN columns shared by all analyses
pecan_cols <- c(
  "pecan_n_edges", "pecan_strength", "pecan_density",
  "DASS_T1_total", "DASS_T2_total"
)

# Helper: add node counts, severity, loops; select final column set.
# main_nar_nodes_col: which column from nc_main to use as nar_nodes for main.
# main_nar_le4_col:   which column from loops_main to use as nar_le4 for main.
build_combined <- function(pilot_nm, main_nm,
                           pilot_nar_cols, main_nar_cols,
                           out_nar_names,
                           main_nar_nodes_col = "nar_nodes",
                           main_nar_le4_col   = "nar_le4") {
  # pilot part: always uses nar_nodes and nar_le4
  pilot_part <- pilot_nm |>
    dplyr::select(vp, study, all_of(pecan_cols), all_of(pilot_nar_cols)) |>
    rename_with(~ out_nar_names, all_of(pilot_nar_cols)) |>
    left_join(nc_pilot |> dplyr::select(vp, nar_nodes, pecan_nodes), by = "vp") |>
    left_join(sev_pilot, by = "vp") |>
    left_join(loops_pilot |> dplyr::select(vp, pecan_le4, pecan_reg_le4, nar_le4),
              by = "vp")

  # main part: pull the correct nar_nodes and nar_le4 variant
  nc_main_sel <- nc_main |>
    transmute(vp, pecan_nodes,
              nar_nodes = .data[[main_nar_nodes_col]])

  loops_main_sel <- loops_main |>
    transmute(vp, pecan_le4, pecan_reg_le4,
              nar_le4 = .data[[main_nar_le4_col]])

  main_part <- main_nm |>
    dplyr::select(vp, study, all_of(pecan_cols), all_of(main_nar_cols)) |>
    rename_with(~ out_nar_names, all_of(main_nar_cols)) |>
    left_join(nc_main_sel,    by = "vp") |>
    left_join(sev_main,       by = "vp") |>
    left_join(loops_main_sel, by = "vp")

  bind_rows(pilot_part, main_part)
}

# -- Analysis 1: pilot nar + main nar_ex --
combined_a1_unreg <- build_combined(
  pilot_nm           = nm_pilot_unreg,
  main_nm            = nm_main_unreg,
  pilot_nar_cols     = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_cols      = c("nar_ex_n_edges", "nar_ex_strength", "nar_ex_density"),
  out_nar_names      = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_nodes_col = "nar_ex_nodes",
  main_nar_le4_col   = "nar_ex_le4"
)

combined_a1_reg <- build_combined(
  pilot_nm           = nm_pilot_reg,
  main_nm            = nm_main_reg,
  pilot_nar_cols     = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_cols      = c("nar_ex_n_edges", "nar_ex_strength", "nar_ex_density"),
  out_nar_names      = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_nodes_col = "nar_ex_nodes",
  main_nar_le4_col   = "nar_ex_le4"
)

# -- Analysis 2: pilot nar + main nar (full: explicit + implicit) --
combined_a2_unreg <- build_combined(
  pilot_nm           = nm_pilot_unreg,
  main_nm            = nm_main_unreg,
  pilot_nar_cols     = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_cols      = c("nar_n_edges", "nar_strength", "nar_density"),
  out_nar_names      = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_nodes_col = "nar_nodes",
  main_nar_le4_col   = "nar_le4"
)

combined_a2_reg <- build_combined(
  pilot_nm           = nm_pilot_reg,
  main_nm            = nm_main_reg,
  pilot_nar_cols     = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_cols      = c("nar_n_edges", "nar_strength", "nar_density"),
  out_nar_names      = c("nar_n_edges", "nar_strength", "nar_density"),
  main_nar_nodes_col = "nar_nodes",
  main_nar_le4_col   = "nar_le4"
)

cat("=== COMBINED DASS-21 EXTERNAL VALIDITY (N =", nrow(combined_a1_unreg), ") ===\n")
cat("Pilot: n =", sum(combined_a1_unreg$study == "pilot"),
    "| Main: n =", sum(combined_a1_unreg$study == "main"),
    "| Total: n =", nrow(combined_a1_unreg), "\n")
cat("T1 DASS — Mean:", round(mean(combined_a1_unreg$DASS_T1_total, na.rm = TRUE), 2),
    "| SD:", round(sd(combined_a1_unreg$DASS_T1_total, na.rm = TRUE), 2),
    "| Range:", min(combined_a1_unreg$DASS_T1_total, na.rm = TRUE), "–",
    max(combined_a1_unreg$DASS_T1_total, na.rm = TRUE), "\n")
cat("T2 DASS — Mean:", round(mean(combined_a1_unreg$DASS_T2_total, na.rm = TRUE), 2),
    "| SD:", round(sd(combined_a1_unreg$DASS_T2_total, na.rm = TRUE), 2),
    "| Missing:", sum(is.na(combined_a1_unreg$DASS_T2_total)), "\n")

# ---------------------------------------------------------------------------
# 2. Spearman correlation helper
# ---------------------------------------------------------------------------

# nar_strength is excluded from both lists: nar edges are binary (all width = 1),
# so nar_strength == nar_n_edges for every participant. Including it would
# duplicate the n_edges correlation under a misleading weighted-metric label.
metrics_unreg <- c(
  "nar_n_edges", "nar_density",
  "nar_nodes", "nar_le4",
  "pecan_n_edges", "pecan_strength", "pecan_density",
  "pecan_nodes", "pecan_mean_severity", "pecan_le4"
)

# Metrics for regulated PECAN (swap pecan_le4 → pecan_reg_le4)
metrics_reg <- c(
  "nar_n_edges", "nar_density",
  "nar_nodes", "nar_le4",
  "pecan_n_edges", "pecan_strength", "pecan_density",
  "pecan_nodes", "pecan_mean_severity", "pecan_reg_le4"
)

ext_cor_combined <- function(data, metrics, dass_col, dass_label) {
  map_dfr(metrics, function(metric) {
    valid <- !is.na(data[[metric]]) & !is.na(data[[dass_col]])
    ct    <- cor.test(data[[metric]][valid], data[[dass_col]][valid],
                      method = "spearman", exact = FALSE)
    tibble(
      dass    = dass_label,
      metric  = metric,
      rho     = round(ct$estimate, 3),
      p_value = round(ct$p.value, 4),
      n       = sum(valid)
    )
  })
}

run_analysis <- function(combined_unreg, combined_reg, label) {
  cat("\n\n===", label, "===\n")

  # Unregulated
  cat("\n--- PECAN unregulated ---\n")
  res_unreg <- bind_rows(
    ext_cor_combined(combined_unreg, metrics_unreg, "DASS_T1_total", "T1"),
    ext_cor_combined(combined_unreg, metrics_unreg, "DASS_T2_total", "T2")
  )
  print(res_unreg, n = Inf)

  # Regulated
  cat("\n--- PECAN regulated ---\n")
  res_reg <- bind_rows(
    ext_cor_combined(combined_reg, metrics_reg, "DASS_T1_total", "T1"),
    ext_cor_combined(combined_reg, metrics_reg, "DASS_T2_total", "T2")
  )
  print(res_reg, n = Inf)

  list(unreg = res_unreg, reg = res_reg,
       data_unreg = combined_unreg, data_reg = combined_reg)
}

# ---------------------------------------------------------------------------
# 3. Run all analyses
# ---------------------------------------------------------------------------

res_a1 <- run_analysis(
  combined_a1_unreg, combined_a1_reg,
  "Analysis 1: pilot nar + main nar_ex (explicit only)"
)

res_a2 <- run_analysis(
  combined_a2_unreg, combined_a2_reg,
  "Analysis 2: pilot nar + main nar (full: explicit + implicit)"
)

# ---------------------------------------------------------------------------
# 4. Save results
# ---------------------------------------------------------------------------

# Analysis 1
write.csv(res_a1$unreg,
          "output/combined/tables/2a_combined_DASS_A1_pilot-nar_main-narex_unreg.csv",
          row.names = FALSE)
write.csv(res_a1$reg,
          "output/combined/tables/2a_combined_DASS_A1_pilot-nar_main-narex_reg.csv",
          row.names = FALSE)
write.csv(res_a1$data_unreg,
          "output/combined/tables/2a_combined_metrics_A1_pilot-nar_main-narex_unreg.csv",
          row.names = FALSE)

# Analysis 2
write.csv(res_a2$unreg,
          "output/combined/tables/2b_combined_DASS_A2_pilot-nar_main-nar_unreg.csv",
          row.names = FALSE)
write.csv(res_a2$reg,
          "output/combined/tables/2b_combined_DASS_A2_pilot-nar_main-nar_reg.csv",
          row.names = FALSE)
write.csv(res_a2$data_unreg,
          "output/combined/tables/2b_combined_metrics_A2_pilot-nar_main-nar_unreg.csv",
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Wide summary table — all analyses side by side
# ---------------------------------------------------------------------------

make_wide <- function(results, analysis, pecan) {
  results |>
    pivot_wider(
      names_from  = dass,
      values_from = c(rho, p_value, n),
      names_glue  = "{dass}_{.value}"
    ) |>
    mutate(analysis = analysis, pecan = pecan, .before = metric)
}

wide_all <- bind_rows(
  make_wide(res_a1$unreg, "A1_pilot-nar+main-narex", "unreg"),
  make_wide(res_a1$reg,   "A1_pilot-nar+main-narex", "reg"),
  make_wide(res_a2$unreg, "A2_pilot-nar+main-nar",   "unreg"),
  make_wide(res_a2$reg,   "A2_pilot-nar+main-nar",   "reg")
)

cat("\n[Summary] All combined DASS correlations — wide format:\n")
print(wide_all, n = Inf)

write.csv(wide_all,
          "output/combined/tables/2c_combined_DASS_all_analyses_wide.csv",
          row.names = FALSE)

cat("\n=== Combined DASS analysis complete. Results saved to output/combined/tables/ ===\n")

# ---------------------------------------------------------------------------
# 5b. BH-corrected external validity correlations
# ---------------------------------------------------------------------------
# Applies Benjamini-Hochberg FDR correction to the Section 3 Spearman
# correlations as a complement to the uncorrected results.
#
# ── TEST FAMILY DEFINITION ──────────────────────────────────────────────────
#
# A test family should contain all tests that jointly address ONE research
# question, so that BH controls the expected proportion of false discoveries
# within that question.
#
# Option A — one family per DASS timepoint [PRIMARY, RECOMMENDED]:
#   Family T1 : all metrics × DASS T1  (10 tests)
#   Family T2 : all metrics × DASS T2  (10 tests)
#   Rationale : T1 (concurrent validity) and T2 (predictive follow-up) are
#   distinct substantive questions. Within each, every metric is a candidate
#   answer to "which network properties associate with distress?", making them
#   a coherent family. nar and PECAN metrics belong together because the same
#   overarching question applies to both instruments.
#
# Option B — one family per instrument × timepoint [ALTERNATIVE]:
#   Family nar-T1   : 4 tests   Family nar-T2   : 4 tests
#   Family PECAN-T1 : 6 tests   Family PECAN-T2 : 6 tests
#   Rationale : nar and PECAN are methodologically distinct; their metrics
#   may be seen as addressing separate sub-questions.
#   Caveat : with only 4 tests per nar family, BH is very lenient there
#   (corrected p ≈ raw p × 4/4 for the most significant test), limiting its
#   practical value. Option A is therefore preferred.
#
# Correction is applied within each analysis variant (A1/A2 × unreg/reg)
# separately — variants are sensitivity checks of each other, not independent
# replications, so they should not be pooled into one family.
# ────────────────────────────────────────────────────────────────────────────

nar_metrics_set <- c("nar_n_edges", "nar_density", "nar_nodes", "nar_le4")

add_bh <- function(df) {
  df |>
    # Option A: family = DASS timepoint
    group_by(dass) |>
    mutate(p_fdr_A = round(p.adjust(p_value, method = "BH"), 4)) |>
    # Option B: family = instrument × timepoint
    group_by(dass, instrument = if_else(metric %in% nar_metrics_set, "nar", "pecan")) |>
    mutate(p_fdr_B = round(p.adjust(p_value, method = "BH"), 4)) |>
    ungroup() |>
    dplyr::select(-instrument)
}

bh_all <- bind_rows(
  res_a1$unreg |> add_bh() |> mutate(analysis = "A1_pilot-nar+main-narex", pecan = "unreg"),
  res_a1$reg   |> add_bh() |> mutate(analysis = "A1_pilot-nar+main-narex", pecan = "reg"),
  res_a2$unreg |> add_bh() |> mutate(analysis = "A2_pilot-nar+main-nar",   pecan = "unreg"),
  res_a2$reg   |> add_bh() |> mutate(analysis = "A2_pilot-nar+main-nar",   pecan = "reg")
) |>
  dplyr::select(analysis, pecan, dass, metric, rho, p_value, p_fdr_A, p_fdr_B, n)

cat("\n=== Section 5b: BH-corrected external validity correlations ===\n")
cat("p_fdr_A : BH within DASS timepoint — 10 tests/family [recommended]\n")
cat("p_fdr_B : BH within instrument x timepoint — 4 (nar) or 6 (PECAN) tests/family\n\n")
print(bh_all, n = Inf)

write.csv(bh_all,
          "output/combined/tables/3_combined_DASS_BH_corrected.csv",
          row.names = FALSE)

cat("\nSaved BH-corrected correlations to combined_DASS_BH_corrected.csv\n")

# ---------------------------------------------------------------------------
# 6. Partial Spearman correlations: all metrics ~ DASS | n_nodes
# ---------------------------------------------------------------------------
# Rationale: network size (n_nodes) may confound any metric that scales with
# the number of nodes (edges, density, loops). Partial Spearman correlations
# isolate each metric's unique association with DASS above and beyond n_nodes.
#
# Method: ppcor::pcor.test(x, y, z, method = "spearman")
#   x = network metric, y = DASS score, z = n_nodes (covariate)
# Each metric is paired with its own node-count covariate:
#   nar metrics  → nar_nodes
#   pecan metrics → pecan_nodes
# ---------------------------------------------------------------------------

library(ppcor)

partial_spearman_fn <- function(data, metric_col, node_col, dass_col,
                                 analysis_label, source_label) {
  ok <- complete.cases(data[[metric_col]], data[[node_col]], data[[dass_col]])
  n  <- sum(ok)
  if (n < 5) {
    return(tibble(analysis = analysis_label, source = source_label,
                  metric = metric_col, node_covariate = node_col,
                  dass = sub("DASS_|_total", "", dass_col),
                  rho_raw = NA_real_, rho_partial = NA_real_,
                  p_value = NA_real_, n = n))
  }
  x <- data[[metric_col]][ok]
  y <- data[[dass_col]][ok]
  z <- data[[node_col]][ok]
  rho_raw <- round(cor(x, y, method = "spearman"), 3)
  pc      <- ppcor::pcor.test(x, y, z, method = "spearman")
  tibble(
    analysis       = analysis_label,
    source         = source_label,
    metric         = metric_col,
    node_covariate = node_col,
    dass           = sub("DASS_|_total", "", dass_col),
    rho_raw        = rho_raw,
    rho_partial    = round(pc$estimate,  3),
    p_value        = round(pc$p.value,   4),
    n              = n
  )
}

# Metric–covariate pairs: nar metrics use nar_nodes; pecan metrics use pecan_nodes
partial_specs_unreg <- list(
  # nar_strength omitted: equals nar_n_edges for binary nar edges
  list(metric = "nar_n_edges",          nodes = "nar_nodes",   source = "nar"),
  list(metric = "nar_density",          nodes = "nar_nodes",   source = "nar"),
  list(metric = "nar_le4",              nodes = "nar_nodes",   source = "nar"),
  list(metric = "pecan_n_edges",        nodes = "pecan_nodes", source = "pecan"),
  list(metric = "pecan_strength",       nodes = "pecan_nodes", source = "pecan"),
  list(metric = "pecan_density",        nodes = "pecan_nodes", source = "pecan"),
  list(metric = "pecan_le4",            nodes = "pecan_nodes", source = "pecan"),
  list(metric = "pecan_mean_severity",  nodes = "pecan_nodes", source = "pecan")
)

# Regulated: swap pecan_le4 → pecan_reg_le4
partial_specs_reg <- partial_specs_unreg
partial_specs_reg[[8]]$metric <- "pecan_reg_le4"

run_partial <- function(data, specs, analysis_label) {
  map_dfr(specs, function(s) {
    map_dfr(c("DASS_T1_total", "DASS_T2_total"), function(dass_col) {
      partial_spearman_fn(data, s$metric, s$nodes, dass_col,
                          analysis_label, s$source)
    })
  })
}

partial_results_all <- bind_rows(
  run_partial(combined_a1_unreg, partial_specs_unreg,
              "A1_pilot-nar+main-narex | PECAN unreg"),
  run_partial(combined_a1_reg,   partial_specs_reg,
              "A1_pilot-nar+main-narex | PECAN reg"),
  run_partial(combined_a2_unreg, partial_specs_unreg,
              "A2_pilot-nar+main-nar | PECAN unreg"),
  run_partial(combined_a2_reg,   partial_specs_reg,
              "A2_pilot-nar+main-nar | PECAN reg")
)

cat("\n=== Partial Spearman: all metrics ~ DASS | n_nodes ===\n")
print(partial_results_all, n = Inf)

write.csv(partial_results_all,
          "output/combined/tables/4_combined_DASS_partial_correlations.csv",
          row.names = FALSE)

cat("\nSaved partial correlation table.\n")

# ---------------------------------------------------------------------------
# 6b. BH-corrected partial Spearman correlations
# ---------------------------------------------------------------------------
# Applies Benjamini-Hochberg FDR correction to the partial correlations from
# Section 6, mirroring the logic of Section 5b.
#
# ── TEST FAMILY DEFINITION ──────────────────────────────────────────────────
#
# Option A — one family per analysis variant × DASS timepoint [RECOMMENDED]:
#   8 metrics per family (3 nar + 5 PECAN; nar_strength excluded as duplicate,
#   node counts excluded as covariates).
#   Rationale: identical to Section 5b — T1 and T2 are distinct questions,
#   and all metrics within a timepoint jointly answer "which network properties
#   uniquely associate with distress beyond network size?"
#
# Option B — one family per instrument × analysis variant × DASS timepoint:
#   Family nar   : 3 tests (n_edges, density, le4)
#   Family PECAN : 5 tests (n_edges, strength, density, le4, mean_severity)
#   Same caveat as Section 5b: nar families are small (3 tests), making BH
#   very lenient. Option A is therefore preferred.
# ────────────────────────────────────────────────────────────────────────────

partial_bh <- partial_results_all |>
  # Option A: family = analysis variant × DASS timepoint
  group_by(analysis, dass) |>
  mutate(p_fdr_A = round(p.adjust(p_value, method = "BH"), 4)) |>
  # Option B: family = instrument × analysis variant × DASS timepoint
  group_by(analysis, dass, source) |>
  mutate(p_fdr_B = round(p.adjust(p_value, method = "BH"), 4)) |>
  ungroup() |>
  dplyr::select(analysis, source, metric, node_covariate, dass,
                rho_raw, rho_partial, p_value, p_fdr_A, p_fdr_B, n)

cat("\n=== Section 6b: BH-corrected partial Spearman correlations ===\n")
cat("p_fdr_A : BH within analysis variant x DASS timepoint — 8 tests/family [recommended]\n")
cat("p_fdr_B : BH within instrument x analysis variant x DASS timepoint — 3 (nar) or 5 (PECAN) tests/family\n\n")
print(partial_bh, n = Inf)

write.csv(partial_bh,
          "output/combined/tables/4b_combined_DASS_partial_BH_corrected.csv",
          row.names = FALSE)

cat("\nSaved BH-corrected partial correlations to combined_DASS_partial_BH_corrected.csv\n")

# ---------------------------------------------------------------------------
# 7. DASS change score: amount of change ~ initial DASS + PECAN node severity
# ---------------------------------------------------------------------------
# DASS_change = T2 - T1. Negative values indicate improvement.
# Research questions:
#   (a) Does the amount of change correlate with initial distress (DASS T1)?
#   (b) Does the amount of change correlate with PECAN mean node severity
#       (i.e., does network severity predict who changes most)?
# Note: pecan_mean_severity is the same across reg/unreg; results will be
# identical for both — both variants are included for completeness.

add_dass_change <- function(df) {
  df |> mutate(DASS_change = DASS_T2_total - DASS_T1_total)
}

combined_a1_unreg_ch <- add_dass_change(combined_a1_unreg)
combined_a1_reg_ch   <- add_dass_change(combined_a1_reg)
combined_a2_unreg_ch <- add_dass_change(combined_a2_unreg)
combined_a2_reg_ch   <- add_dass_change(combined_a2_reg)

cat("\n=== Section 7: DASS change score (T2 - T1) descriptives ===\n")
cat("n with both T1 and T2:", sum(complete.cases(combined_a1_unreg_ch$DASS_change)), "\n")
cat("Mean change:", round(mean(combined_a1_unreg_ch$DASS_change, na.rm = TRUE), 2),
    "| SD:", round(sd(combined_a1_unreg_ch$DASS_change, na.rm = TRUE), 2),
    "| Range:", min(combined_a1_unreg_ch$DASS_change, na.rm = TRUE), "–",
    max(combined_a1_unreg_ch$DASS_change, na.rm = TRUE), "\n")

change_cor_fn <- function(data, analysis_label) {
  vars <- list(
    list(x = "DASS_T1_total",       label = "DASS_T1_initial"),
    list(x = "pecan_mean_severity",  label = "PECAN_mean_node_severity")
  )
  map_dfr(vars, function(v) {
    ok <- complete.cases(data$DASS_change, data[[v$x]])
    if (sum(ok) < 5) {
      return(tibble(analysis = analysis_label, predictor = v$label,
                    rho = NA_real_, p_value = NA_real_, n = sum(ok)))
    }
    ct <- cor.test(data$DASS_change[ok], data[[v$x]][ok],
                   method = "spearman", exact = FALSE)
    tibble(
      analysis  = analysis_label,
      predictor = v$label,
      rho       = round(ct$estimate, 3),
      p_value   = round(ct$p.value,  4),
      n         = sum(ok)
    )
  })
}

dass_change_results <- bind_rows(
  change_cor_fn(combined_a1_unreg_ch, "A1 | PECAN unreg"),
  change_cor_fn(combined_a1_reg_ch,   "A1 | PECAN reg"),
  change_cor_fn(combined_a2_unreg_ch, "A2 | PECAN unreg"),
  change_cor_fn(combined_a2_reg_ch,   "A2 | PECAN reg")
)

cat("\n=== Section 7: Spearman correlations with DASS change (T2 - T1) ===\n")
print(dass_change_results, n = Inf)

write.csv(dass_change_results,
          "output/combined/tables/5_combined_DASS_change_correlations.csv",
          row.names = FALSE)

cat("\nSaved DASS change correlation table.\n")

# ---------------------------------------------------------------------------
# 8. Network density ~ DASS: focused summary
# ---------------------------------------------------------------------------
# Density is already included in the main correlation analysis (Section 3);
# this section extracts and presents the density rows explicitly as a
# standalone result, and saves a dedicated CSV.
#
# Hypothesis: higher network density reflects a more interconnected symptom
# structure and may be associated with higher distress (positive correlation),
# or alternatively, more regulation capacity may lower density while stress
# remains high (negative correlation). Direction is empirical.

density_results <- wide_all |>
  filter(metric %in% c("nar_density", "pecan_density")) |>
  dplyr::select(analysis, pecan, metric,
               T1_rho, T1_p_value, T1_n,
               T2_rho, T2_p_value, T2_n)

cat("\n=== Section 8: Network density ~ DASS correlations ===\n")
print(density_results, n = Inf)

write.csv(density_results,
          "output/combined/tables/6_combined_DASS_density_summary.csv",
          row.names = FALSE)

cat("\nSaved density–DASS summary table.\n")

# ---------------------------------------------------------------------------
# 9a. Steiger's t-test: nar vs. PECAN — do the two methods associate equally
#    strongly with DASS?
# ---------------------------------------------------------------------------
# For each metric with a parallel nar and PECAN version, the Williams (1959) /
# Steiger (1980) t-test for two dependent correlations sharing one variable
# (DASS) is applied:
#
#   r12 = Spearman(nar_metric,   DASS)
#   r13 = Spearman(pecan_metric, DASS)
#   r23 = Spearman(nar_metric,   pecan_metric)  ← inter-method correlation
#
# H0: the nar and PECAN associations with DASS are equal in magnitude.
# H1: they differ.
#
# Returns a t-statistic on n - 3 df and a two-tailed p-value.

# Williams (1959) / Steiger (1980) t-test for two dependent correlations
# sharing one variable (y = DASS). Returns a t-statistic on n-3 df and a
# two-tailed p-value. All values are guaranteed plain numeric scalars.
#
# r12 = cor(nar,   DASS)   r13 = cor(PECAN, DASS)   r23 = cor(nar, PECAN)
steiger_t_test <- function(r12, r13, r23, n) {
  if (any(!is.finite(c(r12, r13, r23)))) return(list(t = NA_real_, p = NA_real_))
  det_R <- 1 - r12^2 - r13^2 - r23^2 + 2 * r12 * r13 * r23
  if (det_R <= 0) return(list(t = NA_real_, p = NA_real_))
  denom <- sqrt(2 * ((n - 1) / (n - 3)) * det_R +
                  ((r12 + r13)^2 / 4) * (1 - r23)^3)
  if (denom == 0) return(list(t = NA_real_, p = NA_real_))
  t_stat <- as.numeric((r12 - r13) * sqrt((n - 1) * (1 + r23)) / denom)
  p_val  <- as.numeric(2 * pt(-abs(t_stat), df = n - 3))
  list(t = t_stat, p = p_val)
}

steiger_fn <- function(data, nar_col, pecan_col, dass_col,
                       analysis_label, dass_label) {
  ok <- complete.cases(data[[nar_col]], data[[pecan_col]], data[[dass_col]])
  n  <- sum(ok)
  if (n < 5) {
    return(tibble(
      analysis = analysis_label, dass = dass_label,
      nar_metric = nar_col, pecan_metric = pecan_col,
      r_nar_dass = NA_real_, r_pecan_dass = NA_real_,
      r_nar_pecan = NA_real_, t_stat = NA_real_, p_value = NA_real_, n = n
    ))
  }
  x  <- data[[nar_col]][ok]
  zz <- data[[pecan_col]][ok]
  y  <- data[[dass_col]][ok]

  r12 <- cor(x,  y,  method = "spearman")   # nar   ~ DASS
  r13 <- cor(zz, y,  method = "spearman")   # PECAN ~ DASS
  r23 <- cor(x,  zz, method = "spearman")   # nar   ~ PECAN (inter-method)

  res <- steiger_t_test(r12, r13, r23, n)

  tibble(
    analysis     = analysis_label,
    dass         = dass_label,
    nar_metric   = nar_col,
    pecan_metric = pecan_col,
    r_nar_dass   = round(r12,    3),
    r_pecan_dass = round(r13,    3),
    r_nar_pecan  = round(r23,    3),
    t_stat       = round(res$t,  3),
    p_value      = round(res$p,  4),
    n            = n
  )
}

run_steiger <- function(data, analysis_label, pairs) {
  map_dfr(pairs, function(pr) {
    map_dfr(c("DASS_T1_total", "DASS_T2_total"), function(dc) {
      steiger_fn(data, pr$nar, pr$pecan, dc,
                 analysis_label, sub("DASS_|_total", "", dc))
    })
  })
}

# Parallel metric pairs: (nar column, PECAN column).
# Loops use pecan_le4 for unregulated and pecan_reg_le4 for regulated analyses.
steiger_pairs_unreg <- list(
  list(nar = "nar_n_edges", pecan = "pecan_n_edges"),
  list(nar = "nar_density",  pecan = "pecan_density"),
  list(nar = "nar_nodes",    pecan = "pecan_nodes"),
  list(nar = "nar_le4",      pecan = "pecan_le4")
)

steiger_pairs_reg <- list(
  list(nar = "nar_n_edges", pecan = "pecan_n_edges"),
  list(nar = "nar_density",  pecan = "pecan_density"),
  list(nar = "nar_nodes",    pecan = "pecan_nodes"),
  list(nar = "nar_le4",      pecan = "pecan_reg_le4")
)

# FDR correction applied within each analysis × DASS-timepoint cell (across
# the 4 metric pairs). Bonferroni would be inappropriate: the 4 analysis
# variants are sensitivity checks, and the metric pairs are correlated
# (density is a function of n_edges and n_nodes).
steiger_results <- bind_rows(
  run_steiger(combined_a1_unreg, "A1 | PECAN unreg", steiger_pairs_unreg),
  run_steiger(combined_a1_reg,   "A1 | PECAN reg",   steiger_pairs_reg),
  run_steiger(combined_a2_unreg, "A2 | PECAN unreg", steiger_pairs_unreg),
  run_steiger(combined_a2_reg,   "A2 | PECAN reg",   steiger_pairs_reg)
) |>
  group_by(analysis, dass) |>
  mutate(p_fdr = round(p.adjust(p_value, method = "BH"), 4)) |>
  ungroup()

cat("\n=== Section 9: Steiger's t-test — nar vs. PECAN DASS associations ===\n")
cat("H0: nar and PECAN correlate equally strongly with DASS (df = n - 3).\n\n")
print(steiger_results, n = Inf)

write.csv(steiger_results,
          "output/combined/tables/7_combined_DASS_steiger_nar_vs_pecan.csv",
          row.names = FALSE)

cat("\nSaved Steiger's t-test results (with FDR correction) to",
    "combined_DASS_steiger_nar_vs_pecan.csv\n")

# ---------------------------------------------------------------------------
# 9b. Steiger's t-test: PECAN unregulated vs. regulated
# ---------------------------------------------------------------------------
# Tests whether the regulated and unregulated PECAN networks associate
# equally strongly with DASS. Regulation changes edge counts and density
# but leaves the node set unchanged, so node count is excluded (r23 = 1
# would make the test undefined).
# Source: combined_a1_unreg (PECAN unreg) joined with combined_a1_reg (PECAN
# reg). A1 and A2 carry identical PECAN data, so A1 is used throughout.

pecan_reg_unreg_data <- combined_a1_unreg |>
  dplyr::select(vp, DASS_T1_total, DASS_T2_total,
                pecan_n_edges_unreg = pecan_n_edges,
                pecan_density_unreg = pecan_density,
                pecan_le4_unreg     = pecan_le4) |>
  left_join(
    combined_a1_reg |>
      dplyr::select(vp,
                    pecan_n_edges_reg = pecan_n_edges,
                    pecan_density_reg = pecan_density,
                    pecan_le4_reg     = pecan_reg_le4),
    by = "vp"
  )

steiger_pairs_pecan_reg <- list(
  list(nar = "pecan_n_edges_unreg", pecan = "pecan_n_edges_reg"),
  list(nar = "pecan_density_unreg",  pecan = "pecan_density_reg"),
  list(nar = "pecan_le4_unreg",      pecan = "pecan_le4_reg")
)

steiger_pecan_reg_unreg <- map_dfr(steiger_pairs_pecan_reg, function(pr) {
  map_dfr(c("DASS_T1_total", "DASS_T2_total"), function(dc) {
    steiger_fn(pecan_reg_unreg_data, pr$nar, pr$pecan, dc,
               "PECAN unreg vs reg", sub("DASS_|_total", "", dc))
  })
}) |>
  group_by(dass) |>
  mutate(p_fdr = round(p.adjust(p_value, method = "BH"), 4)) |>
  ungroup()

cat("\n=== Section 9b: Steiger's t-test — PECAN unreg vs. PECAN reg ===\n")
cat("(r_nar_dass = PECAN unreg ~ DASS; r_pecan_dass = PECAN reg ~ DASS)\n\n")
print(steiger_pecan_reg_unreg, n = Inf)

write.csv(steiger_pecan_reg_unreg,
          "output/combined/tables/7b_combined_DASS_steiger_pecan_unreg_vs_reg.csv",
          row.names = FALSE)
cat("\nSaved PECAN unreg vs. reg Steiger results.\n")

# ---------------------------------------------------------------------------
# 9c. Steiger's t-test: nar (full) vs. nar_ex (explicit only)
# ---------------------------------------------------------------------------
# Applicable only to main study participants — the pilot study has no
# explicit/implicit split.
# Tests whether including implicit connections (nar) vs. restricting to
# explicit connections only (nar_ex) changes the association with DASS.
# Built directly from nm_main_unreg which carries both column sets.

nar_vs_ex_data <- nm_main_unreg |>
  dplyr::select(vp, DASS_T1_total, DASS_T2_total,
                nar_n_edges, nar_density,
                nar_ex_n_edges, nar_ex_density) |>
  left_join(nc_main    |> dplyr::select(vp, nar_nodes, nar_ex_nodes), by = "vp") |>
  left_join(loops_main |> dplyr::select(vp, nar_le4,  nar_ex_le4),   by = "vp")

steiger_pairs_narex <- list(
  list(nar = "nar_n_edges", pecan = "nar_ex_n_edges"),
  list(nar = "nar_density",  pecan = "nar_ex_density"),
  list(nar = "nar_nodes",    pecan = "nar_ex_nodes"),
  list(nar = "nar_le4",      pecan = "nar_ex_le4")
)

steiger_nar_vs_ex <- map_dfr(steiger_pairs_narex, function(pr) {
  map_dfr(c("DASS_T1_total", "DASS_T2_total"), function(dc) {
    steiger_fn(nar_vs_ex_data, pr$nar, pr$pecan, dc,
               "nar vs nar_ex (main only)", sub("DASS_|_total", "", dc))
  })
}) |>
  group_by(dass) |>
  mutate(p_fdr = round(p.adjust(p_value, method = "BH"), 4)) |>
  ungroup()

cat("\n=== Section 9c: Steiger's t-test — nar vs. nar_ex (main study only) ===\n")
cat("(r_nar_dass = nar ~ DASS; r_pecan_dass = nar_ex ~ DASS)\n\n")
cat("Note: nar_nodes vs nar_ex_nodes returns NA if node sets are identical",
    "(r23 = 1 makes the test undefined).\n\n")
print(steiger_nar_vs_ex, n = Inf)

write.csv(steiger_nar_vs_ex,
          "output/combined/tables/7c_combined_DASS_steiger_nar_vs_narex.csv",
          row.names = FALSE)
cat("\nSaved nar vs. nar_ex Steiger results.\n")

# ---------------------------------------------------------------------------
# 10. Node count ~ density correlation (artefact check)
# ---------------------------------------------------------------------------
# Density = n_edges / (n_nodes × (n_nodes − 1)), so density is mathematically
# constrained by node count. If node count is positively associated with DASS
# and negatively associated with density, the density ~ DASS association may
# be a mathematical artefact rather than a genuine structural signal.
#
# This section checks Spearman(n_nodes, density) for each network type.
# Interpret alongside the main DASS correlations:
#   - strong negative rho(nodes, density) + positive rho(nodes, DASS)
#     → the negative density ~ DASS association is likely artefactual.
#   - weak or absent rho(nodes, density)
#     → density carries information beyond node count.

node_density_fn <- function(data, nodes_col, density_col, label) {
  ok <- complete.cases(data[[nodes_col]], data[[density_col]])
  ct <- cor.test(data[[nodes_col]][ok], data[[density_col]][ok],
                 method = "spearman", exact = FALSE)
  tibble(
    network     = label,
    nodes_col   = nodes_col,
    density_col = density_col,    rho         = round(ct$estimate, 3),
    p_value     = round(ct$p.value,  4),
    n           = sum(ok)
  )
}

node_density_results <- bind_rows(
  # nar: A1 (pilot nar + main nar_ex) and A2 (pilot nar + main nar full)
  node_density_fn(combined_a1_unreg, "nar_nodes",   "nar_density",   "nar — A1 (pilot nar + main nar_ex)"),
  node_density_fn(combined_a2_unreg, "nar_nodes",   "nar_density",   "nar — A2 (pilot nar + main nar)"),
  # PECAN: identical across A1/A2; unreg and reg differ in edge set
  node_density_fn(combined_a1_unreg, "pecan_nodes", "pecan_density", "PECAN unregulated"),
  node_density_fn(combined_a1_reg,   "pecan_nodes", "pecan_density", "PECAN regulated")
)

cat("\n=== Section 10: Node count ~ density correlation (artefact check) ===\n")
cat("A strong negative rho, combined with a positive node count ~ DASS\n")
cat("association (see main results), would indicate that the negative\n")
cat("density ~ DASS association is at least partly a mathematical artefact.\n\n")
print(node_density_results)

write.csv(node_density_results,
          "output/combined/tables/8_combined_node_density_artefact_check.csv",
          row.names = FALSE)
cat("\nSaved node count ~ density artefact check to",
    "combined_node_density_artefact_check.csv\n")

# ---------------------------------------------------------------------------
# 10b. Node count ~ DASS correlation (completing the artefact chain)
# ---------------------------------------------------------------------------
# If n_nodes is positively associated with DASS and negatively with density,
# then density ~ DASS is at least partly a mathematical artefact.
# This section checks Spearman(n_nodes, DASS_T1) and Spearman(n_nodes, DASS_T2).

node_dass_fn <- function(data, nodes_col, dass_col, label) {
  ok <- complete.cases(data[[nodes_col]], data[[dass_col]])
  ct <- cor.test(data[[nodes_col]][ok], data[[dass_col]][ok],
                 method = "spearman", exact = FALSE)
  tibble(
    network   = label,
    nodes_col = nodes_col,
    dass      = sub("DASS_|_total", "", dass_col),
    rho       = round(ct$estimate, 3),
    p_value   = round(ct$p.value,  4),
    n         = sum(ok)
  )
}

node_dass_results <- bind_rows(
  node_dass_fn(combined_a1_unreg, "nar_nodes",   "DASS_T1_total", "nar — A1"),
  node_dass_fn(combined_a1_unreg, "nar_nodes",   "DASS_T2_total", "nar — A1"),
  node_dass_fn(combined_a2_unreg, "nar_nodes",   "DASS_T2_total", "nar — A2"),
  node_dass_fn(combined_a2_unreg, "nar_nodes",   "DASS_T1_total", "nar — A2"),
  node_dass_fn(combined_a1_unreg, "pecan_nodes", "DASS_T1_total", "PECAN"),
  node_dass_fn(combined_a1_unreg, "pecan_nodes", "DASS_T2_total", "PECAN")
)

cat("\n=== Section 10b: Node count ~ DASS correlation ===\n")
cat("A positive rho here, combined with the negative rho(nodes, density)\n")
cat("from Section 10, completes the artefact chain for density ~ DASS.\n\n")
print(node_dass_results)

write.csv(node_dass_results,
          "output/combined/tables/8b_combined_node_dass_artefact_check.csv",
          row.names = FALSE)
cat("\nSaved node count ~ DASS artefact check to",
    "combined_node_dass_artefact_check.csv\n")

# ---------------------------------------------------------------------------
# 11. DASS change analysis
# ---------------------------------------------------------------------------
# DASS dropped substantially from T1 to T2. This raises two questions:
#   (a) Is the improvement predicted by T1 level? (regression to the mean)
#   (b) Do network metrics predict the amount of change, independently of T1?
#
# DASS_change = DASS_T1_total - DASS_T2_total  (positive = improvement)
#
# 11a: Descriptives of change scores + Spearman(DASS_T1, DASS_change)
# 11b: Partial Spearman of each metric with DASS_change, controlling for
#      n_nodes (same as main analysis) and additionally for DASS_T1
#      (i.e., do metrics predict change BEYOND what T1 already explains?)

# Helper: partial Spearman with two covariates (nodes + T1)
partial_spearman_change_fn <- function(data, metric_col, node_col,
                                       analysis_label, source_label) {
  data <- data |>
    mutate(DASS_change = DASS_T1_total - DASS_T2_total)

  ok <- complete.cases(data[[metric_col]], data[[node_col]],
                       data$DASS_change, data$DASS_T1_total)
  n  <- sum(ok)
  if (n < 6) {
    return(tibble(
      analysis = analysis_label, source = source_label,
      metric = metric_col,
      rho_raw = NA_real_, rho_partial_nodes = NA_real_,
      rho_partial_nodes_t1 = NA_real_,
      p_nodes = NA_real_, p_nodes_t1 = NA_real_, n = n
    ))
  }
  x        <- data[[metric_col]][ok]
  y        <- data$DASS_change[ok]
  z_nodes  <- data[[node_col]][ok]
  z_t1     <- data$DASS_T1_total[ok]

  rho_raw    <- round(cor(x, y, method = "spearman"), 3)
  pc_nodes   <- ppcor::pcor.test(x, y, z_nodes,                    method = "spearman")
  pc_both    <- ppcor::pcor.test(x, y, data.frame(z_nodes, z_t1),  method = "spearman")

  tibble(
    analysis             = analysis_label,
    source               = source_label,
    metric               = metric_col,
    rho_raw              = rho_raw,
    rho_partial_nodes    = round(pc_nodes$estimate, 3),
    p_nodes              = round(pc_nodes$p.value,  4),
    rho_partial_nodes_t1 = round(pc_both$estimate,  3),
    p_nodes_t1           = round(pc_both$p.value,   4),
    n                    = n
  )
}

run_change_analysis <- function(data, analysis_label, pecan_nodes_col,
                                nar_nodes_col, nar_metrics, pecan_metrics) {
  d <- data |> mutate(DASS_change = DASS_T1_total - DASS_T2_total)

  # 11a: change score descriptives + T1 ~ change
  ok <- complete.cases(d$DASS_change, d$DASS_T1_total)
  ct <- cor.test(d$DASS_T1_total[ok], d$DASS_change[ok],
                 method = "spearman", exact = FALSE)
  cat(sprintf(
    "\n--- %s ---\n  DASS change (T1-T2): Mean = %.2f, SD = %.2f, Range = %.0f – %.0f (n = %d)\n",
    analysis_label,
    mean(d$DASS_change, na.rm = TRUE), sd(d$DASS_change, na.rm = TRUE),
    min(d$DASS_change, na.rm = TRUE),  max(d$DASS_change, na.rm = TRUE),
    sum(!is.na(d$DASS_change))
  ))
  cat(sprintf("  Spearman(DASS_T1, DASS_change): rho = %.3f, p = %.4f (n = %d)\n",
              ct$estimate, ct$p.value, sum(ok)))

  # 11b: metric ~ change partial correlations
  bind_rows(
    map_dfr(nar_metrics,   ~ partial_spearman_change_fn(data, .x, nar_nodes_col,
                                                         analysis_label, "nar")),
    map_dfr(pecan_metrics, ~ partial_spearman_change_fn(data, .x, pecan_nodes_col,
                                                         analysis_label, "pecan"))
  )
}

# Key metrics: all main partial-correlation metrics, applied to change score
change_nar_metrics   <- c("nar_n_edges", "nar_density", "nar_le4")
change_pecan_metrics <- c("pecan_n_edges", "pecan_strength", "pecan_density",
                          "pecan_le4", "pecan_mean_severity")

cat("\n=== Section 11: DASS change analysis ===\n")
cat("DASS_change = T1 - T2 (positive = improvement).\n")
cat("Partial Spearman: rho_partial_nodes controls for n_nodes only;\n")
cat("rho_partial_nodes_t1 additionally controls for DASS_T1\n")
cat("(tests whether metrics predict change beyond baseline severity).\n")

change_results <- bind_rows(
  run_change_analysis(combined_a1_unreg, "A1 | PECAN unreg",
                      "pecan_nodes", "nar_nodes",
                      change_nar_metrics, change_pecan_metrics),
  run_change_analysis(combined_a1_reg,   "A1 | PECAN reg",
                      "pecan_nodes", "nar_nodes",
                      change_nar_metrics, change_pecan_metrics),
  run_change_analysis(combined_a2_unreg, "A2 | PECAN unreg",
                      "pecan_nodes", "nar_nodes",
                      change_nar_metrics, change_pecan_metrics),
  run_change_analysis(combined_a2_reg,   "A2 | PECAN reg",
                      "pecan_nodes", "nar_nodes",
                      change_nar_metrics, change_pecan_metrics)
)

cat("\n--- Full results table ---\n")
print(change_results, n = Inf)

write.csv(change_results,
          "output/combined/tables/9_combined_DASS_change_analysis.csv",
          row.names = FALSE)
cat("\nSaved DASS change analysis to combined_DASS_change_analysis.csv\n")
