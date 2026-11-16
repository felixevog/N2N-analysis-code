# =============================================================================
# N2N Pilot Study — Master's Thesis Data Analysis
# Comparing Narrative-to-Network (N2N/nar) with PECAN (gold standard)
# N = 22 participants (VP01–VP22)
# =============================================================================

# ===========================================================================
# SECTION 0: SETUP & DATA LOADING
# ===========================================================================

rm(list = ls())   # clear environment so all objects are rebuilt from source files

library(tidyverse)
library(psych)
library(qgraph)
library(PECAN2)
library(irr)
library(ppcor)
suppressWarnings(library(afex))

# Restore default warning behaviour (warn = 1: print warnings, do not error).
# This prevents warn = 2 sessions from converting warnings to hard errors.
options(warn = 1)

# --- 0.1 Create working directories and copy raw data to data/bearbeitet/ ---
# IMPORTANT: All analysis reads from data/bearbeitet/. Raw files are never modified.

dir.create("data/bearbeitet/nar",          recursive = TRUE, showWarnings = FALSE)
dir.create("data/bearbeitet/PECAN",        recursive = TRUE, showWarnings = FALSE)
dir.create("output/pilotstudy/tables",     recursive = TRUE, showWarnings = FALSE)
dir.create("output/pilotstudy/plots",      recursive = TRUE, showWarnings = FALSE)
dir.create("output/pilotstudy/networks",   recursive = TRUE, showWarnings = FALSE)

file.copy(
  "data/rawdata/pilotstudy/N2N T1 Data pilotstudy.csv",
  "data/bearbeitet/N2N T1 Data pilotstudy.csv",
  overwrite = FALSE
)

nar_raw_pilot   <- list.files("data/rawdata/pilotstudy/nar",   pattern = "\\.csv$", full.names = TRUE)
pecan_raw_pilot <- list.files("data/rawdata/pilotstudy/PECAN", pattern = "\\.csv$", full.names = TRUE)
file.copy(nar_raw_pilot,   "data/bearbeitet/nar/",   overwrite = TRUE)
file.copy(pecan_raw_pilot, "data/bearbeitet/PECAN/", overwrite = TRUE)
file.copy(
  "data/rawdata/pilotstudy/N2N T2 Data pilotstudy.csv",
  "data/bearbeitet/N2N T2 Data pilotstudy.csv",
  overwrite = FALSE
)

# --- 0.2 Load T1 survey data ---
# File is UTF-16LE encoded and tab-separated (exported from SoSciSurvey)

t1_pilot <- read.delim(
  "data/bearbeitet/N2N T1 Data pilotstudy.csv",
  fileEncoding = "UTF-16LE",
  sep          = "\t"
)

# --- 0.3 Load per-participant network data into named lists ---

vps_pilot <- sprintf("VP%02d", 1:22)

nar_nodes_list_pilot   <- setNames(
  lapply(vps_pilot, \(vp) read.csv(paste0("data/bearbeitet/nar/",   vp, "_nar_nodes.csv"))),
  vps_pilot
)
nar_edges_list_pilot   <- setNames(
  lapply(vps_pilot, \(vp) read.csv(paste0("data/bearbeitet/nar/",   vp, "_nar_edges.csv"))),
  vps_pilot
)
pecan_nodes_list_pilot <- setNames(
  lapply(vps_pilot, \(vp) read.csv(paste0("data/bearbeitet/PECAN/", vp, "_PECAN_nodes.csv"))),
  vps_pilot
)
pecan_edges_list_pilot <- setNames(
  lapply(vps_pilot, \(vp) read.csv(paste0("data/bearbeitet/PECAN/", vp, "_PECAN_edges.csv"))),
  vps_pilot
)

# r2 data: one folder per rater, auto-detected from rawdata/pilotstudy/nar_r2.*/
# Each folder is named e.g. nar_r2.1_G and contains *_nar_nodes_r2.csv / *_nar_edges_r2.csv.
# r2_raters_pilot: named list (folder name → list(vps_pilot, nodes, edges)); empty rater folders are skipped.
r2_rater_dirs_pilot <- list.dirs(
  "data/rawdata/pilotstudy",
  recursive = FALSE,
  full.names = FALSE
)

r2_rater_dirs_pilot <- sort(
  r2_rater_dirs_pilot[grepl("^nar_r2\\.", r2_rater_dirs_pilot)]
)

r2_raters_pilot <- setNames(lapply(r2_rater_dirs_pilot, function(rd) {
  folder <- file.path("data/rawdata/pilotstudy", rd)
  
  vps_pilot_in_folder <- list.files(
    folder,
    pattern = "_nar_nodes_r2\\.csv$"
  ) |>
    sub("_nar_nodes_r2\\.csv$", "", x = _)
  
  if (length(vps_pilot_in_folder) == 0) return(NULL)
  
  list(
    vps_pilot = vps_pilot_in_folder,
    
    nodes = setNames(
      lapply(vps_pilot_in_folder, \(vp) {
        read.csv(file.path(folder, paste0(vp, "_nar_nodes_r2.csv")))
      }),
      vps_pilot_in_folder
    ),
    
    edges = setNames(
      lapply(vps_pilot_in_folder, \(vp) {
        read.csv(file.path(folder, paste0(vp, "_nar_edges_r2.csv")))
      }),
      vps_pilot_in_folder
    )
  )
}), r2_rater_dirs_pilot)

r2_raters_pilot <- Filter(Negate(is.null), r2_raters_pilot)

cat("r2 rater folders found (n =", length(r2_raters_pilot), "):\n")

for (rd in names(r2_raters_pilot)) {
  cat(
    "  ", rd, ": n =", length(r2_raters_pilot[[rd]]$vps_pilot), "VPs —",
    paste(r2_raters_pilot[[rd]]$vps_pilot, collapse = ", "),
    "\n"
  )
}

n_vps_pilot <- length(vps_pilot)

# Pool all nodes and edges across participants (used in Sections 3b-ALT and 5)
all_nar_nodes_pilot   <- map_dfr(vps_pilot, ~ nar_nodes_list_pilot[[.x]]   |> mutate(vp = .x))
all_nar_edges_pilot   <- map_dfr(vps_pilot, ~ nar_edges_list_pilot[[.x]]   |> mutate(vp = .x))
all_pecan_nodes_pilot <- map_dfr(vps_pilot, ~ pecan_nodes_list_pilot[[.x]] |> mutate(vp = .x, label = str_squish(label)))
all_pecan_edges_pilot <- map_dfr(vps_pilot, ~ pecan_edges_list_pilot[[.x]] |> mutate(vp = .x))

# --- 0.5 PECAN edge regulation ---
# Regulation rule (matching pecanVis edges_regulation):
#   reg_by    = "width"  → rank edges by width, strongest first
#   reg_type  = "nodes"  → cut-off is relative to number of nodes
#   reg_value = 1.0      → keep ≤ n_nodes edges per participant
# Ties at the cut-off boundary are always kept in full.

regulate_pecan_edges <- function(edges_df, nodes_df, reg_value = 1.0) {
  n_nodes <- nrow(nodes_df)
  n_keep  <- floor(n_nodes * reg_value)
  if (nrow(edges_df) == 0 || n_keep <= 0) return(edges_df[integer(0), ])
  edges_sorted <- edges_df |> arrange(desc(width))
  if (nrow(edges_sorted) <= n_keep) return(edges_sorted)
  cutoff_width <- edges_sorted$width[n_keep]
  edges_sorted |> filter(width >= cutoff_width)
}

pecan_edges_reg_list_pilot <- setNames(
  lapply(vps_pilot, \(vp) regulate_pecan_edges(pecan_edges_list_pilot[[vp]], pecan_nodes_list_pilot[[vp]])),
  vps_pilot
)
all_pecan_edges_reg_pilot <- map_dfr(vps_pilot, ~ pecan_edges_reg_list_pilot[[.x]] |> mutate(vp = .x))

# --- 0.4 Helper functions ---

# Jaccard similarity between two sets
jaccard <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  uni   <- length(union(set_a, set_b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# Build a square binary adjacency matrix from an edge list over a given node set --> for the Edge Concordance
#differences build_adj and make_adj: build_adj is binary (0 = edge not present, 1 = edge present), make_adj takes edge weight into account
#build_adj is used in IRR section (comparing two raters' edge patterns via Cohen's kappa), where you only need presence/absence. 
#make_adj is used in section 3b for the Spearman edge concordance, where the same binary logic applies but the weighted option exists for the sensitivity checks.
#output from make_adj ist adj_nar_list_pilot and adj_pecan_list_pilot 
#You'd need to run at minimum through to line ~560 for build_adj_list_pilot to appear, and through to ~845 for adj_nar_list_pilot.
make_adj <- function(edges_df, node_ids, weight_col = NULL) {
  node_ids <- as.character(node_ids)
  n        <- length(node_ids)
  mat      <- matrix(0, n, n, dimnames = list(node_ids, node_ids))
  
  edges_df <- edges_df |>
    mutate(from = as.character(from), to = as.character(to)) |>
    filter(from %in% node_ids, to %in% node_ids) |>
    # Remove duplicate edges (can occur in nar files)
    distinct(from, to, .keep_all = TRUE)
  
  if (nrow(edges_df) == 0) return(mat)
  
  for (i in seq_len(nrow(edges_df))) {
    f  <- edges_df$from[i]
    t_ <- edges_df$to[i]
    mat[f, t_] <- if (!is.null(weight_col)) edges_df[[weight_col]][i] else 1
  }
  mat
}

# --- 0.6 Additional helper functions ---

# Catch errors from any IRR function (covers warnings-as-errors when
# options(warn = 2) is set, as well as genuine errors)
safe_irr <- function(expr) tryCatch(expr, error = function(e) NA_real_)

# Compute in- and out-degree over a given set of node IDs
#counts how many in/outdegree connections a node has for the nodes that are shared between two raters
#so it is used for IRR, specifically for comparing Outdegree Centrality Concordance between two raters

deg_fn <- function(edges_df, ids) {
  out <- edges_df |>
    mutate(from = as.character(from)) |>
    filter(from %in% ids) |>
    count(from, name = "out_deg") |>
    rename(id = from)
  inn <- edges_df |>
    mutate(to = as.character(to)) |>
    filter(to %in% ids) |>
    count(to, name = "in_deg") |>
    rename(id = to)
  tibble(id = ids) |>
    left_join(out, by = "id") |>
    left_join(inn, by = "id") |>
    mutate(across(c(out_deg, in_deg), ~ replace_na(.x, 0L)))
}

# Returns list(icc, alpha); NA if no variation on either side
icc_alpha_fn <- function(v1, v2) {
  if (length(v1) < 2 || sd(v1) == 0 || sd(v2) == 0) return(list(icc = NA_real_, alpha = NA_real_))
  mat <- data.frame(r1 = v1, r2 = v2)
  icc_val   <- safe_irr(round(irr::icc(mat, model = "twoway", type = "consistency")$value, 3))
  alpha_val <- safe_irr(
    withCallingHandlers(
      round(psych::alpha(mat, warnings = FALSE)$total$raw_alpha, 3),
      warning = function(w) invokeRestart("muffleWarning")
    )
  )
  list(icc = icc_val, alpha = alpha_val)
}

# Format mean ± SD with n valid for console summaries
fmt <- function(x) paste0(round(mean(x, na.rm = TRUE), 3),
                           " \u00b1 ", round(sd(x, na.rm = TRUE), 3),
                           "  (n valid = ", sum(!is.na(x)), ")")

# Append Mean / SD / Min / Max summary rows to a tibble before writing CSV.
# For combined files (with a group_col), summaries are computed per group.
append_summary <- function(df, rho_cols, group_col = NULL) {
  stats <- list(Mean = mean, SD = sd, Min = min, Max = max)
  smry  <- map_dfr(names(stats), function(stat) {
    fn <- stats[[stat]]
    base <- if (!is.null(group_col)) group_by(df, across(all_of(group_col))) else df
    base |>
      summarise(across(all_of(rho_cols), ~ round(fn(.x, na.rm = TRUE), 3)),
                .groups = "drop") |>
      mutate(vp = paste0("--- ", stat, " ---"))
  })
  bind_rows(df, smry)
}

# Spearman rho + p for two numeric vectors (per-participant centrality)
rho_fn <- function(x, y) {
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  list(rho = round(ct$estimate, 3), p = round(ct$p.value, 4))
}

# Format range as "min – max"
rng <- function(x) paste0(round(min(x, na.rm = TRUE), 3), " \u2013 ", round(max(x, na.rm = TRUE), 3))

# Pre-compute co-occurrence: for each VP, which directed node pairs co-occur?
#all possible edges for each VP are listed --> used to compare the actual edges later on
build_cooccurrence <- function(nodes_list, vps_pilot) {
  map_dfr(vps_pilot, function(vp) {
    ids <- nodes_list[[vp]]$id
    if (length(ids) < 2) return(tibble())
    expand.grid(from = ids, to = ids, stringsAsFactors = FALSE) |>
      as_tibble() |>
      filter(from != to) |>
      mutate(vp = vp)
  })
}

# Matrix trace helper for cycle counting
#helper for identifying and counting feedback loops for each VP
#used to create nar_cycle_list_pilot and pecan_cycle_list_pilot
mat_trace <- function(M) sum(diag(M))

# Count directed cycles of length 2, 3, 4 via matrix powers
count_cycles_le4 <- function(adj) {
  A  <- adj
  A2 <- A  %*% A
  A3 <- A2 %*% A
  A4 <- A2 %*% A2
  d2 <- diag(A2)

  n2 <- round(mat_trace(A2) / 2)
  n3 <- round(mat_trace(A3) / 3)
  n4 <- round((mat_trace(A4) - 2 * sum(d2^2) + mat_trace(A2)) / 4)

  list(n2 = n2, n3 = n3, n4 = n4, total_le4 = n2 + n3 + n4)
}

# Build binary adjacency matrix from an edge list for a given node set
#differences build_adj and make_adj: build_adj is binary (0 = edge not present, 1 = edge present), make_adj takes edge weight into account
#build_adj is used in IRR section (comparing two raters' edge patterns via Cohen's kappa), where you only need presence/absence. 
#make_adj is used in section 3b for the Spearman edge concordance, where the same binary logic applies but the weighted option exists for the sensitivity checks.
#You'd need to run at minimum through to line ~560 for build_adj_list_pilot to appear, and through to ~845 for adj_nar_list_pilot.
build_adj <- function(edges_df, node_ids) {
  node_ids <- as.character(node_ids)
  n        <- length(node_ids)
  adj      <- matrix(0L, n, n, dimnames = list(node_ids, node_ids))
  edges_clean <- edges_df |>
    mutate(from = as.character(from), to = as.character(to)) |>
    filter(from %in% node_ids, to %in% node_ids) |>
    distinct(from, to)
  for (i in seq_len(nrow(edges_clean))) {
    adj[edges_clean$from[i], edges_clean$to[i]] <- 1L
  }
  adj
}

# Canonical cycle label: lexicographically smallest cyclic rotation
canonical_cycle <- function(labels_vec) {
  n    <- length(labels_vec)
  rots <- sapply(seq_len(n), function(i)
    paste(labels_vec[c(i:n, seq_len(i - 1L))], collapse = " \u2192 "))
  rots[order(rots)[1L]]
}

# Enumerate all directed cycles of length 2–4 with node labels
enumerate_cycles_labeled <- function(edges_df, nodes_df) {
  lbl <- setNames(as.character(nodes_df$label), as.character(nodes_df$id))
  ids <- names(lbl)
  n   <- length(ids)
  if (n < 2) return(character(0))

  adj <- matrix(0L, n, n, dimnames = list(ids, ids))
  edges_clean <- edges_df |>
    mutate(from = as.character(from), to = as.character(to)) |>
    filter(from %in% ids, to %in% ids) |>
    distinct(from, to)
  for (i in seq_len(nrow(edges_clean)))
    adj[edges_clean$from[i], edges_clean$to[i]] <- 1L

  cycles <- character(0)

  # Length-2: mutual pairs
  for (i in seq_len(n))
    for (j in seq_len(n))
      if (i < j && adj[i, j] == 1L && adj[j, i] == 1L)
        cycles <- c(cycles, canonical_cycle(lbl[c(ids[i], ids[j])]))

  # Length-3: i→j→k→i (all distinct)
  for (i in seq_len(n))
    for (j in seq_len(n)) { if (j == i || adj[i,j] == 0L) next
      for (k in seq_len(n)) { if (k == i || k == j || adj[j,k] == 0L || adj[k,i] == 0L) next
        cycles <- c(cycles, canonical_cycle(lbl[c(ids[i], ids[j], ids[k])]))
      }
    }

  # Length-4: i→j→k→l→i (all distinct)
  for (i in seq_len(n))
    for (j in seq_len(n)) { if (j == i || adj[i,j] == 0L) next
      for (k in seq_len(n)) { if (k == i || k == j || adj[j,k] == 0L) next
        for (l in seq_len(n)) { if (l == i || l == j || l == k || adj[k,l] == 0L || adj[l,i] == 0L) next
          cycles <- c(cycles, canonical_cycle(lbl[c(ids[i], ids[j], ids[k], ids[l])]))
        }
      }
    }

  unique(cycles)
}


# ===========================================================================
# SECTION 1: DESCRIPTIVE SAMPLE DATA
# ===========================================================================
#relevant for the master thesis were the following output tables:
#1_descriptives_age.csv
#1_descriptives_gender.csv
#1_descriptives_education.csv

cat("\n=== SECTION 1: DESCRIPTIVE SAMPLE DATA ===\n")

# Age (SD02_01)
age_stats_pilot <- t1_pilot |>
  summarise(
    mean_age = round(mean(SD02_01, na.rm = TRUE), 2),
    sd_age   = round(sd(SD02_01,   na.rm = TRUE), 2),
    min_age  = min(SD02_01, na.rm = TRUE),
    max_age  = max(SD02_01, na.rm = TRUE)
  )
cat("Age — Mean:", age_stats_pilot$mean_age, "| SD:", age_stats_pilot$sd_age,
    "| Range:", age_stats_pilot$min_age, "–", age_stats_pilot$max_age, "\n")

# Gender (SD05): 1 = female, 2 = male
gender_table_pilot <- t1_pilot |>
  mutate(gender = case_when(
    SD05 == 1 ~ "female",
    SD05 == 2 ~ "male",
    TRUE      ~ "other/missing"
  )) |>
  count(gender) |>
  mutate(pct = round(n / sum(n) * 100, 1))
cat("\nGender distribution:\n")
print(gender_table_pilot)

# Education (SD11)
education_table_pilot <- t1_pilot |>
  count(SD11) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  arrange(SD11)
cat("\nEducation distribution (SD11):\n")
print(education_table_pilot)

write.csv(as.data.frame(age_stats_pilot), "output/pilotstudy/tables/1_descriptives_age.csv",       row.names = FALSE)
write.csv(gender_table_pilot,             "output/pilotstudy/tables/1_descriptives_gender.csv",    row.names = FALSE)
write.csv(education_table_pilot,          "output/pilotstudy/tables/1_descriptives_education.csv", row.names = FALSE)


# ===========================================================================
# SECTION 2: RELIABILITY & CONCORDANCE
# ===========================================================================
#relevant for the master thesis were the following output tables:
#2b_cronbach_alpha_DASS21 #internal consistency of the DASS-21
#2c_irr_global_kappa_nar_r2.1_G #interrater-reliability measured between main rater and second rater 1 with a reference set pooled across all participants
#2c_irr_global_kappa_nar_r2.2_A #interrater-reliability measured between main rater and second rater 2 with a reference set pooled across all participants

cat("\n=== SECTION 2: RELIABILITY & CONCORDANCE ===\n")

# --- 2a. DASS-21 sum score (recode 1–4 → 0–3 per standard scoring) ---

dass_items_pilot <- paste0("DA01_", sprintf("%02d", 1:21))

dass_data_pilot <- t1_pilot |>
  dplyr::select(VP01_01, all_of(dass_items_pilot)) |>
  # Standardise VP ID capitalisation (raw data has "Vp15" instead of "VP15")
  mutate(VP01_01 = toupper(VP01_01)) |>
  # Recode: original 1–4 → standard 0–3
  mutate(across(all_of(dass_items_pilot), ~ .x - 1)) |>
  mutate(DASS_total = rowSums(across(all_of(dass_items_pilot)), na.rm = FALSE))

cat("DASS-21 sum scores per participant:\n")
print(dass_data_pilot |> dplyr::select(VP01_01, DASS_total))
cat("DASS-21 total — Mean:", round(mean(dass_data_pilot$DASS_total, na.rm = TRUE), 2),
    "| SD:", round(sd(dass_data_pilot$DASS_total, na.rm = TRUE), 2),
    "| Range:", min(dass_data_pilot$DASS_total, na.rm = TRUE), "–",
    max(dass_data_pilot$DASS_total, na.rm = TRUE), "\n")

write.csv(dass_data_pilot |> dplyr::select(VP01_01, DASS_total),
          "output/pilotstudy/tables/2a_DASS21_scores.csv", row.names = FALSE)

# --- 2a-T2. DASS-21 T2 sum scores ---
# T2 file was copied to bearbeitet/ in Section 0.1

t2_pilot <- read.delim(
  "data/bearbeitet/N2N T2 Data pilotstudy.csv",
  fileEncoding = "UTF-16LE",
  sep = "\t"
)

dass_t2_items_pilot <- paste0("DA01_", sprintf("%02d", 1:21))
dass_t2_items_data_pilot <- t2_pilot |>
  dplyr::select(vp = VP, all_of(dass_t2_items_pilot)) |>
  mutate(across(all_of(dass_t2_items_pilot), ~ .x - 1))

dass_t2_data_pilot  <- dass_t2_items_data_pilot |>
  mutate(DASS_T2_total = rowSums(across(all_of(dass_t2_items_pilot)), na.rm = FALSE)) |>
  dplyr::select(vp, DASS_T2_total)

cat("\nDASS-21 T2 sum scores:\n")
print(dass_t2_data_pilot)
cat("T2 DASS — Mean:", round(mean(dass_t2_data_pilot$DASS_T2_total, na.rm = TRUE), 2),
    "| SD:", round(sd(dass_t2_data_pilot$DASS_T2_total, na.rm = TRUE), 2),
    "| Range:", min(dass_t2_data_pilot$DASS_T2_total, na.rm = TRUE), "–",
    max(dass_t2_data_pilot$DASS_T2_total, na.rm = TRUE),
    "| Missing:", sum(is.na(dass_t2_data_pilot$DASS_T2_total)), "\n")

write.csv(dass_t2_data_pilot, "output/pilotstudy/tables/2a_DASS21_T2_scores.csv", row.names = FALSE)

# --- 2a-OV. DASS-21 descriptives overview (T1 and T2) ---

dass_overview_pilot <- bind_rows(
  tibble(
    timepoint = "T1",
    n         = sum(!is.na(dass_data_pilot$DASS_total)),
    M         = round(mean(dass_data_pilot$DASS_total,      na.rm = TRUE), 2),
    SD        = round(sd(dass_data_pilot$DASS_total,        na.rm = TRUE), 2),
    Min       = min(dass_data_pilot$DASS_total,             na.rm = TRUE),
    Max       = max(dass_data_pilot$DASS_total,             na.rm = TRUE)
  ),
  tibble(
    timepoint = "T2",
    n         = sum(!is.na(dass_t2_data_pilot$DASS_T2_total)),
    M         = round(mean(dass_t2_data_pilot$DASS_T2_total, na.rm = TRUE), 2),
    SD        = round(sd(dass_t2_data_pilot$DASS_T2_total,   na.rm = TRUE), 2),
    Min       = min(dass_t2_data_pilot$DASS_T2_total,        na.rm = TRUE),
    Max       = max(dass_t2_data_pilot$DASS_T2_total,        na.rm = TRUE)
  )
)

cat("\nDASS-21 total score overview (scale 0–63):\n")
print(dass_overview_pilot)

write.csv(dass_overview_pilot,
          "output/pilotstudy/tables/2a_DASS21_descriptives_overview.csv",
          row.names = FALSE)

# --- 2b. Cronbach's alpha for DASS-21 (T1 and T2) ---

alpha_result_pilot <- withCallingHandlers(
  psych::alpha(dass_data_pilot |> dplyr::select(all_of(dass_items_pilot)), check.keys = TRUE),
  warning = function(w) invokeRestart("muffleWarning")
)

alpha_result_t2_pilot <- withCallingHandlers(
  psych::alpha(dass_t2_items_data_pilot |> dplyr::select(all_of(dass_t2_items_pilot)), check.keys = TRUE),
  warning = function(w) invokeRestart("muffleWarning")
)

cat("\nCronbach's Alpha for DASS-21 — T1:\n")
print(alpha_result_pilot$total)
cat("\nCronbach's Alpha for DASS-21 — T2:\n")
print(alpha_result_t2_pilot$total)

alpha_combined_pilot <- bind_rows(
  as.data.frame(alpha_result_pilot$total)    |> mutate(timepoint = "T1", .before = 1),
  as.data.frame(alpha_result_t2_pilot$total) |> mutate(timepoint = "T2", .before = 1)
)
write.csv(alpha_combined_pilot, "output/pilotstudy/tables/2b_cronbach_alpha_DASS21.csv",
          row.names = FALSE)

# --- 2c. Inter-rater Reliability (nar r1 vs. nar r2) — looped per rater ---
#
# r2_raters_pilot is a named list built in Section 0 from rawdata/pilotstudy/nar_r2.*/ folders.
# Each rater is processed independently so that IRR estimates are not conflated across raters.
#
# Per rater, two sets of metrics are computed:
#   (i)  Nodes & edges (binary presence/absence) → Cohen's Kappa (irr::kappa2)
#         + Jaccard index (node sets)
#         + ICC two-way consistency & Cronbach's alpha for in-/out-degree
#   (ii) Kappa with a within-rater global reference set
#         (node/edge vocab = union across all VPs of that rater)
#
# Per-rater CSVs: 2c_irr_<rater_id>.csv, 2c_irr_global_kappa_<rater_id>.csv
# Combined CSVs:  2c_irr_all_raters.csv, 2c_irr_global_kappa_all_raters.csv

if (length(r2_raters_pilot) == 0) {
  cat("\nSection 2c skipped: no r2 rater data available for pilot study.\n")
} else {

all_irr_pilot        <- list()
all_global_pilot     <- list()
deg_fn_list_pilot    <- list()
build_adj_list_pilot <- list()
global_ref_list_pilot <- list()

for (rater_id in names(r2_raters_pilot)) {

  vps_r2            <- r2_raters_pilot[[rater_id]]$vps_pilot
  nar_nodes_r2_list_pilot <- r2_raters_pilot[[rater_id]]$nodes
  nar_edges_r2_list_pilot <- r2_raters_pilot[[rater_id]]$edges

  cat(sprintf("\n=== IRR: Rater %s (n = %d VPs: %s) ===\n",
              rater_id, length(vps_r2), paste(vps_r2, collapse = ", ")))

  # --- (i) Main: kappa + Jaccard + ICC/alpha ---
  irr_results_pilot <- map_dfr(vps_r2, function(vp) {
    n1 <- nar_nodes_list_pilot[[vp]]
    e1 <- distinct(nar_edges_list_pilot[[vp]], from, to)
    n2 <- nar_nodes_r2_list_pilot[[vp]]
    e2 <- distinct(nar_edges_r2_list_pilot[[vp]], from, to)

    ids1       <- as.character(n1$id)
    ids2       <- as.character(n2$id)
    all_ids    <- sort(union(ids1, ids2))
    shared_ids <- intersect(ids1, ids2)

    # Nodes: Cohen's Kappa on binary presence vectors (union vocab)
    node_r1 <- as.integer(all_ids %in% ids1)
    node_r2 <- as.integer(all_ids %in% ids2)
    kappa_nodes <- if (sd(node_r1) == 0 && sd(node_r2) == 0) NA_real_ else
      safe_irr(round(irr::kappa2(data.frame(r1 = node_r1, r2 = node_r2))$value, 3))

    # Edges: Cohen's Kappa on flattened binary adjacency vectors
    adj1    <- build_adj(e1, all_ids)
    adj2    <- build_adj(e2, all_ids)
    idx     <- which(row(adj1) != col(adj1))
    edge_r1 <- adj1[idx]; edge_r2 <- adj2[idx]
    kappa_edges <- if (sd(edge_r1) == 0 && sd(edge_r2) == 0) NA_real_ else
      safe_irr(round(irr::kappa2(data.frame(r1 = edge_r1, r2 = edge_r2))$value, 3))

    # Degree: ICC + Cronbach's alpha over shared nodes
    if (length(shared_ids) >= 2) {
      deg1    <- deg_fn(e1, shared_ids)
      deg2    <- deg_fn(e2, shared_ids)
      res_in  <- icc_alpha_fn(deg1$in_deg,  deg2$in_deg)
      res_out <- icc_alpha_fn(deg1$out_deg, deg2$out_deg)
    } else {
      res_in  <- list(icc = NA_real_, alpha = NA_real_)
      res_out <- list(icc = NA_real_, alpha = NA_real_)
    }

    tibble(
      vp             = vp,
      n_nodes_r1     = length(ids1),
      n_nodes_r2     = length(ids2),
      n_shared_nodes = length(shared_ids),
      kappa_nodes    = kappa_nodes,
      kappa_edges    = kappa_edges,
      jaccard_nodes  = jaccard(ids1, ids2),
      icc_indeg      = res_in$icc,
      alpha_indeg    = res_in$alpha,
      icc_outdeg     = res_out$icc,
      alpha_outdeg   = res_out$alpha
    )
  })

  # Save deg_fn output per VP for inspection
  # Save build_adj output per VP for inspection
  build_adj_list_pilot[[rater_id]] <- setNames(
    lapply(vps_r2, function(vp) {
      e1      <- distinct(nar_edges_list_pilot[[vp]], from, to)
      e2      <- distinct(nar_edges_r2_list_pilot[[vp]], from, to)
      all_ids <- sort(union(
        as.character(nar_nodes_list_pilot[[vp]]$id),
        as.character(nar_nodes_r2_list_pilot[[vp]]$id)
      ))
      list(
        rater1 = build_adj(e1, all_ids),
        rater2 = build_adj(e2, all_ids)
      )
    }),
    vps_r2
  )

  deg_fn_list_pilot[[rater_id]] <- setNames(
    lapply(vps_r2, function(vp) {
      e1         <- distinct(nar_edges_list_pilot[[vp]], from, to)
      e2         <- distinct(nar_edges_r2_list_pilot[[vp]], from, to)
      shared_ids <- intersect(
        as.character(nar_nodes_list_pilot[[vp]]$id),
        as.character(nar_nodes_r2_list_pilot[[vp]]$id)
      )
      if (length(shared_ids) < 2) return(NULL)
      list(
        rater1 = deg_fn(e1, shared_ids),
        rater2 = deg_fn(e2, shared_ids)
      )
    }),
    vps_r2
  )

  cat("\nIRR (kappa / ICC / Jaccard) —", rater_id, ":\n")
  print(irr_results_pilot)
  cat("\n--- Summary (mean ± SD | Range) ---\n")
  cat("Nodes  — Kappa:       ", fmt(irr_results_pilot$kappa_nodes),   "| Range:", rng(irr_results_pilot$kappa_nodes),   "\n")
  cat("Edges  — Kappa:       ", fmt(irr_results_pilot$kappa_edges),   "| Range:", rng(irr_results_pilot$kappa_edges),   "\n")
  cat("Nodes  — Jaccard:     ", fmt(irr_results_pilot$jaccard_nodes), "| Range:", rng(irr_results_pilot$jaccard_nodes), "\n")
  cat("In-deg — ICC:         ", fmt(irr_results_pilot$icc_indeg),     "| Range:", rng(irr_results_pilot$icc_indeg),     "\n")
  cat("In-deg — Cronbach α:  ", fmt(irr_results_pilot$alpha_indeg),   "| Range:", rng(irr_results_pilot$alpha_indeg),   "\n")
  cat("Out-deg — ICC:        ", fmt(irr_results_pilot$icc_outdeg),    "| Range:", rng(irr_results_pilot$icc_outdeg),    "\n")
  cat("Out-deg — Cronbach α: ", fmt(irr_results_pilot$alpha_outdeg),  "| Range:", rng(irr_results_pilot$alpha_outdeg),  "\n")
  write.csv(
    append_summary(irr_results_pilot,
      c("kappa_nodes", "kappa_edges", "jaccard_nodes",
        "icc_indeg", "alpha_indeg", "icc_outdeg", "alpha_outdeg")),
    sprintf("output/pilotstudy/tables/2c_irr_%s.csv", rater_id), row.names = FALSE
  )

  # --- (ii) Kappa with within-rater global reference set ---
  # Global vocab = union across all VPs rated by this rater.
  # Makes kappa values structurally comparable within the rater's VP set.
  global_node_ids_r_pilot <- sort(unique(c(
    unlist(lapply(vps_r2, \(vp) as.character(nar_nodes_list_pilot[[vp]]$id))),
    unlist(lapply(vps_r2, \(vp) as.character(nar_nodes_r2_list_pilot[[vp]]$id)))
  )))
  global_edge_ids_r_pilot <- sort(unique(c(
    unlist(lapply(vps_r2, \(vp) {
      e <- distinct(nar_edges_list_pilot[[vp]], from, to)
      paste(e$from, e$to, sep = "->")
    })),
    unlist(lapply(vps_r2, \(vp) {
      e <- distinct(nar_edges_r2_list_pilot[[vp]], from, to)
      paste(e$from, e$to, sep = "->")
    }))
  )))

  # Save reference set for inspection
  global_ref_list_pilot[[rater_id]] <- list(
    node_ids = global_node_ids_r_pilot,
    edge_ids = global_edge_ids_r_pilot,
    # Binary presence vectors per VP (what actually goes into kappa)
    vectors  = setNames(
      lapply(vps_r2, function(vp) {
        ids1 <- as.character(nar_nodes_list_pilot[[vp]]$id)
        ids2 <- as.character(nar_nodes_r2_list_pilot[[vp]]$id)
        e1   <- distinct(nar_edges_list_pilot[[vp]], from, to)
        e2   <- distinct(nar_edges_r2_list_pilot[[vp]], from, to)
        ep1  <- paste(e1$from, e1$to, sep = "->")
        ep2  <- paste(e2$from, e2$to, sep = "->")
        list(
          gnode_r1 = as.integer(global_node_ids_r_pilot %in% ids1),
          gnode_r2 = as.integer(global_node_ids_r_pilot %in% ids2),
          gedge_r1 = as.integer(global_edge_ids_r_pilot %in% ep1),
          gedge_r2 = as.integer(global_edge_ids_r_pilot %in% ep2)
        )
      }),
      vps_r2
    )
  )

  irr_global_kappa_pilot <- map_dfr(vps_r2, function(vp) {
    ids1 <- as.character(nar_nodes_list_pilot[[vp]]$id)
    ids2 <- as.character(nar_nodes_r2_list_pilot[[vp]]$id)
    e1   <- distinct(nar_edges_list_pilot[[vp]], from, to)
    e2   <- distinct(nar_edges_r2_list_pilot[[vp]], from, to)
    ep1  <- paste(e1$from, e1$to, sep = "->")
    ep2  <- paste(e2$from, e2$to, sep = "->")

    gnode_r1 <- as.integer(global_node_ids_r_pilot %in% ids1)
    gnode_r2 <- as.integer(global_node_ids_r_pilot %in% ids2)
    kappa_nodes_global <- if (sd(gnode_r1) == 0 && sd(gnode_r2) == 0) NA_real_ else
      safe_irr(round(irr::kappa2(data.frame(r1 = gnode_r1, r2 = gnode_r2))$value, 3))

    gedge_r1 <- as.integer(global_edge_ids_r_pilot %in% ep1)
    gedge_r2 <- as.integer(global_edge_ids_r_pilot %in% ep2)
    kappa_edges_global <- if (sd(gedge_r1) == 0 && sd(gedge_r2) == 0) NA_real_ else
      safe_irr(round(irr::kappa2(data.frame(r1 = gedge_r1, r2 = gedge_r2))$value, 3))

    tibble(vp = vp, kappa_nodes_global = kappa_nodes_global, kappa_edges_global = kappa_edges_global)
  })

  cat("\nIRR (global ref kappa) —", rater_id, ":\n")
  print(irr_global_kappa_pilot)
  cat("Nodes — Kappa (global ref):", fmt(irr_global_kappa_pilot$kappa_nodes_global),
      "| Range:", rng(irr_global_kappa_pilot$kappa_nodes_global), "\n")
  cat("Edges — Kappa (global ref):", fmt(irr_global_kappa_pilot$kappa_edges_global),
      "| Range:", rng(irr_global_kappa_pilot$kappa_edges_global), "\n")
  cat("(Reference: n =", length(global_node_ids_r_pilot), "node IDs,",
      length(global_edge_ids_r_pilot), "directed edge pairs)\n")
  write.csv(
    append_summary(irr_global_kappa_pilot, c("kappa_nodes_global", "kappa_edges_global")),
    sprintf("output/pilotstudy/tables/2c_irr_global_kappa_%s.csv", rater_id), row.names = FALSE
  )

  # Collect for cross-rater combined tables
  all_irr_pilot[[rater_id]]    <- irr_results_pilot      |> mutate(rater = rater_id, .before = 1)
  all_global_pilot[[rater_id]] <- irr_global_kappa_pilot |> mutate(rater = rater_id, .before = 1)

} # end for rater_id

# --- Cross-rater combined tables ---
irr_all_raters_pilot <- bind_rows(all_irr_pilot)
irr_all_global_pilot <- bind_rows(all_global_pilot)

cat("\n=== Combined IRR: mean per rater (kappa / ICC / Jaccard) ===\n")
irr_all_raters_pilot |>
  group_by(rater) |>
  summarise(
    n_vps_pilot       = n(),
    kappa_nodes = round(mean(kappa_nodes,   na.rm = TRUE), 3),
    kappa_edges = round(mean(kappa_edges,   na.rm = TRUE), 3),
    jaccard     = round(mean(jaccard_nodes, na.rm = TRUE), 3),
    icc_indeg   = round(mean(icc_indeg,     na.rm = TRUE), 3),
    icc_outdeg  = round(mean(icc_outdeg,    na.rm = TRUE), 3),
    .groups = "drop"
  ) |> print()

cat("\n=== Combined IRR: overall mean, SD, and range across all raters ===\n")
irr_summary_metrics_pilot <- c("kappa_nodes", "kappa_edges", "jaccard_nodes",
                          "icc_indeg",   "alpha_indeg",
                          "icc_outdeg",  "alpha_outdeg")
irr_overall_summary_pilot <- map_dfr(irr_summary_metrics_pilot, function(m) {
  x <- irr_all_raters_pilot[[m]][!is.na(irr_all_raters_pilot[[m]])]
  tibble(metric = m,
         n      = length(x),
         mean   = round(mean(x), 3),
         sd     = round(sd(x),   3),
         min    = round(min(x),  3),
         max    = round(max(x),  3))
})
print(irr_overall_summary_pilot)

write.csv(irr_overall_summary_pilot,
  "output/pilotstudy/tables/2c_irr_overall_summary.csv", row.names = FALSE)
write.csv(irr_all_raters_pilot,
  "output/pilotstudy/tables/2c_irr_all_raters.csv", row.names = FALSE)
write.csv(irr_all_global_pilot,
  "output/pilotstudy/tables/2c_irr_global_kappa_all_raters.csv", row.names = FALSE)

cat("\n=== Combined IRR: overall mean, SD, and range — global reference kappa ===\n")
irr_global_overall_summary_pilot <- map_dfr(
  c("kappa_nodes_global", "kappa_edges_global"),
  function(m) {
    x <- irr_all_global_pilot[[m]][!is.na(irr_all_global_pilot[[m]])]
    tibble(metric = m,
           n      = length(x),
           mean   = round(mean(x), 3),
           sd     = round(sd(x),   3),
           min    = round(min(x),  3),
           max    = round(max(x),  3))
  }
)
print(irr_global_overall_summary_pilot)
write.csv(irr_global_overall_summary_pilot,
  "output/pilotstudy/tables/2c_irr_global_kappa_overall_summary.csv", row.names = FALSE)

} # end if (length(r2_raters_pilot) > 0)

# ===========================================================================
# SECTION 3: STRUCTURAL VALIDITY
# ===========================================================================
#relevant for the master thesis were the following output tables:
#3a_jaccard_summary #jaccard indices of node overlap between nar network and PECAN network
#3b_edge_concordance_spearman_combined #summary of Spearman correlations between narrative network and PECAN networks
#3c_centrality_spearman_combined #summary of centrality Spearman correlations between narrative network and PECAN networks
#3c_top_node_agreement_overview #asks the question whether nar and PECAN networks identified the same nodes to have the highest in/outdegree

cat("\n=== SECTION 3: STRUCTURAL VALIDITY ===\n")

# --- 3a. Jaccard Index (node concordance per participant) ---

jaccard_results_pilot <- tibble(
  vp      = vps_pilot,
  n_nar   = map_int(vps_pilot, ~ nrow(nar_nodes_list_pilot[[.x]])),
  n_pecan = map_int(vps_pilot, ~ nrow(pecan_nodes_list_pilot[[.x]])),
  n_shared = map_int(vps_pilot, ~ length(intersect(
    nar_nodes_list_pilot[[.x]]$id, pecan_nodes_list_pilot[[.x]]$id
  ))),
  jaccard = map_dbl(vps_pilot, ~ jaccard(
    nar_nodes_list_pilot[[.x]]$id, pecan_nodes_list_pilot[[.x]]$id
  ))
) |>
  mutate(
    # Proportion of nar nodes that are also in PECAN (nar coverage by PECAN)
    prop_shared_of_nar   = round(n_shared / n_nar,   3),
    # Proportion of PECAN nodes that are also in nar (PECAN coverage by nar)
    prop_shared_of_pecan = round(n_shared / n_pecan, 3)
  )

jaccard_summary_pilot <- jaccard_results_pilot |>
  summarise(
    mean_J              = round(mean(jaccard,              na.rm = TRUE), 3),
    sd_J                = round(sd(jaccard,                na.rm = TRUE), 3),
    min_J               = round(min(jaccard,               na.rm = TRUE), 3),
    max_J               = round(max(jaccard,               na.rm = TRUE), 3),
    mean_n_shared       = round(mean(n_shared,             na.rm = TRUE), 3),
    sd_n_shared         = round(sd(n_shared,               na.rm = TRUE), 3),
    min_n_shared        = round(min(n_shared,              na.rm = TRUE), 3),
    max_n_shared        = round(max(n_shared,              na.rm = TRUE), 3),
    mean_prop_of_nar    = round(mean(prop_shared_of_nar,   na.rm = TRUE), 3),
    mean_prop_of_pecan  = round(mean(prop_shared_of_pecan, na.rm = TRUE), 3)
  )

cat("\nJaccard Index per participant (nar vs. PECAN node concordance):\n")
print(jaccard_results_pilot)
cat("\nJaccard summary:\n")
print(jaccard_summary_pilot)

write.csv(jaccard_results_pilot, "output/pilotstudy/tables/3a_jaccard_per_participant.csv", row.names = FALSE)
write.csv(jaccard_summary_pilot, "output/pilotstudy/tables/3a_jaccard_summary.csv", row.names = FALSE)


# --- 3b. Edge concordance (Spearman on binary adjacency matrices) ---
# For each VP: build binary adjacency matrices over concordant (shared) nodes,
# flatten to vectors, compute Spearman rho.

# Pre-build adjacency matrix lists for inspection / reuse
adj_nar_list_pilot <- setNames(
  lapply(vps_pilot, function(vp) {
    shared_ids <- intersect(nar_nodes_list_pilot[[vp]]$id, pecan_nodes_list_pilot[[vp]]$id)
    if (length(shared_ids) < 2) return(NULL)
    make_adj(nar_edges_list_pilot[[vp]], shared_ids)
  }),
  vps_pilot
)

adj_pecan_list_pilot <- setNames(
  lapply(vps_pilot, function(vp) {
    shared_ids <- intersect(nar_nodes_list_pilot[[vp]]$id, pecan_nodes_list_pilot[[vp]]$id)
    if (length(shared_ids) < 2) return(NULL)
    mat <- make_adj(pecan_edges_list_pilot[[vp]], shared_ids)
    mat[mat > 0] <- 1
    mat
  }),
  vps_pilot
)

# Save flattened edge vectors (off-diagonal of adjacency matrices) per VP for inspection
edge_vec_list_pilot <- setNames(
  lapply(vps_pilot, function(vp) {
    if (is.null(adj_nar_list_pilot[[vp]])) return(NULL)
    idx <- which(row(adj_nar_list_pilot[[vp]]) != col(adj_nar_list_pilot[[vp]]))
    list(
      nar_vec   = adj_nar_list_pilot[[vp]][idx],
      pecan_vec = adj_pecan_list_pilot[[vp]][idx]
    )
  }),
  vps_pilot
)

edge_cor_results_pilot <- map_dfr(vps_pilot, function(vp) {
  shared_ids <- intersect(nar_nodes_list_pilot[[vp]]$id, pecan_nodes_list_pilot[[vp]]$id)
  
  if (length(shared_ids) < 3) {
    return(tibble(vp = vp, rho = NA_real_, p = NA_real_, n_shared_nodes = length(shared_ids)))
  }
  
  adj_nar   <- make_adj(nar_edges_list_pilot[[vp]],   shared_ids)
 
   # Dichotomize PECAN edges to 0/1
  adj_pecan <- make_adj(pecan_edges_list_pilot[[vp]], shared_ids)
  adj_pecan[adj_pecan > 0] <- 1
  
#Flatten all off-diagonal cells to a vector (directed: both upper and lower triangle)
#Hier werden alle Zellen außerhalb der Hauptdiagonale herausgenommen und in Vektoren umgewandelt. Die Diagonale wird ausgeschlossen, weil sie Selbstkanten wie A → A darstellen würde.
#Bei drei Knoten gibt es beispielsweise sechs mögliche gerichtete Kanten:
#A → B, A → C,
#B → A, B → C,
#C → A, C → B
#Die beiden Vektoren könnten so aussehen:
#NAR:   1, 0, 0, 1, 0, 0
#PECAN: 1, 0, 1, 0, 0, 0
  idx <- which(row(adj_nar) != col(adj_nar))
  nar_vec   <- adj_nar[idx]
  pecan_vec <- adj_pecan[idx]
  
  if (sd(nar_vec) == 0 || sd(pecan_vec) == 0) {
    return(tibble(vp = vp, rho = NA_real_, p = NA_real_, n_shared_nodes = length(shared_ids)))
  }
  
  ct <- cor.test(nar_vec, pecan_vec, method = "spearman", exact = FALSE)
  tibble(vp = vp, rho = round(ct$estimate, 3), p = round(ct$p.value, 4),
         n_shared_nodes = length(shared_ids))
})

cat("\nEdge concordance (Spearman rho, binary adjacency matrices):\n")
print(edge_cor_results_pilot)
cat("Mean rho:", round(mean(edge_cor_results_pilot$rho, na.rm = TRUE), 3),
    "| SD:", round(sd(edge_cor_results_pilot$rho,   na.rm = TRUE), 3),
    "| Range:", round(min(edge_cor_results_pilot$rho, na.rm = TRUE), 3),
    "–", round(max(edge_cor_results_pilot$rho, na.rm = TRUE), 3), "\n")

# One-sample Wilcoxon: is mean rho significantly > 0?
cat("\nOne-sample Wilcoxon test (H0: rho = 0, H1: rho > 0):\n")
print(wilcox.test(edge_cor_results_pilot$rho, mu = 0, alternative = "greater", exact = FALSE))
write.csv(
  append_summary(edge_cor_results_pilot, "rho"),
  "output/pilotstudy/tables/3b_edge_concordance_spearman_unreg.csv", row.names = FALSE
)

# --- 3b (REGULATED). Edge concordance on regulated PECAN networks ---

edge_cor_results_reg_pilot <- map_dfr(vps_pilot, function(vp) {
  shared_ids <- intersect(nar_nodes_list_pilot[[vp]]$id, pecan_nodes_list_pilot[[vp]]$id)

  if (length(shared_ids) < 3) {
    return(tibble(vp = vp, rho = NA_real_, p = NA_real_, n_shared_nodes = length(shared_ids)))
  }

  adj_nar   <- make_adj(nar_edges_list_pilot[[vp]], shared_ids)
  adj_pecan <- make_adj(pecan_edges_reg_list_pilot[[vp]], shared_ids)
  adj_pecan[adj_pecan > 0] <- 1

  idx       <- which(row(adj_nar) != col(adj_nar))
  nar_vec   <- adj_nar[idx]
  pecan_vec <- adj_pecan[idx]

  if (sd(nar_vec) == 0 || sd(pecan_vec) == 0) {
    return(tibble(vp = vp, rho = NA_real_, p = NA_real_, n_shared_nodes = length(shared_ids)))
  }

  ct <- cor.test(nar_vec, pecan_vec, method = "spearman", exact = FALSE)
  tibble(vp = vp, rho = round(ct$estimate, 3), p = round(ct$p.value, 4),
         n_shared_nodes = length(shared_ids))
})

cat("\n[REGULATED] Edge concordance (Spearman rho, binary adjacency matrices):\n")
print(edge_cor_results_reg_pilot)
cat("Mean rho:", round(mean(edge_cor_results_reg_pilot$rho, na.rm = TRUE), 3),
    "| SD:", round(sd(edge_cor_results_reg_pilot$rho,   na.rm = TRUE), 3),
    "| Range:", round(min(edge_cor_results_reg_pilot$rho, na.rm = TRUE), 3),
    "–", round(max(edge_cor_results_reg_pilot$rho, na.rm = TRUE), 3), "\n")

cat("\n[REGULATED] One-sample Wilcoxon test (H0: rho = 0, H1: rho > 0):\n")
print(wilcox.test(edge_cor_results_reg_pilot$rho, mu = 0, alternative = "greater", exact = FALSE))

write.csv(
  append_summary(edge_cor_results_reg_pilot, "rho"),
  "output/pilotstudy/tables/3b_edge_concordance_spearman_reg.csv", row.names = FALSE
)

# Combined unregulated + regulated edge concordance
edge_cor_combined_pilot <- bind_rows(
  edge_cor_results_pilot     |> mutate(pecan_version = "unregulated"),
  edge_cor_results_reg_pilot |> mutate(pecan_version = "regulated")
) |> dplyr::select(pecan_version, everything())
write.csv(
  append_summary(edge_cor_combined_pilot, "rho", group_col = "pecan_version"),
  "output/pilotstudy/tables/3b_edge_concordance_spearman_combined.csv", row.names = FALSE
)

# --- 3c. Node centrality comparison (Spearman) ---
# Primary comparison: nar simple degree (in+out, unweighted) vs. PECAN unweighted degree
# (in_connect + out_connect). This is a symmetrical structural comparison.
#
# Secondary: also compute against PECAN weighted degree (sum of widths) and
# sevwod_combined (severity-weighted), reported as sensitivity checks.
# NOTE: Only participants with ≥ 3 shared nodes yield valid correlations;
# with n = 3–4 nodes per correlation, estimates are highly unstable.

nar_degree_list_pilot <- list()

centrality_cor_results_pilot <- map_dfr(vps_pilot, function(vp) {
  pn <- pecan_nodes_list_pilot[[vp]]
  pe <- pecan_edges_list_pilot[[vp]]
  ne <- nar_edges_list_pilot[[vp]]

  na_row <- tibble(vp = vp,
                   rho_in_unwt  = NA_real_, p_in_unwt  = NA_real_, rho_in_wt  = NA_real_,
                   rho_out_unwt = NA_real_, p_out_unwt = NA_real_, rho_out_wt = NA_real_,
                   rho_sevwod   = NA_real_, n_shared   = NA_integer_)

  pecan_cen <- tryCatch(
    pecanCen(pe, pn, centrality_by = "width", sev_weighted = "value"),
    error = function(e) { warning(paste("pecanCen failed for", vp, ":", e$message)); NULL }
  )
  if (is.null(pecan_cen)) return(na_row)

  pecan_cen <- pecan_cen |> mutate(id = as.character(id))

  # nar in-degree and out-degree kept separate
  #nar_degree_list_pilot and nar_degree_reg_list_pilot looking at the centrality for each node/vp to check 
  nar_out    <- ne |> count(from, name = "out_deg") |> rename(id = from)
  nar_in     <- ne |> count(to,   name = "in_deg")  |> rename(id = to)
  nar_degree <- full_join(nar_out, nar_in, by = "id") |>
    mutate(across(c(out_deg, in_deg), ~ replace_na(.x, 0L)),
           id = as.character(id))

  nar_degree_list_pilot[[vp]] <<- nar_degree

  shared <- inner_join(
    nar_degree |> dplyr::select(id, in_deg, out_deg),
    pecan_cen  |> dplyr::select(id, in_connect, out_connect, in_degree, out_degree, sevwod_combined),
    by = "id"
  )

  if (nrow(shared) < 3) return(na_row |> mutate(n_shared = nrow(shared)))

  r_in_unwt  <- rho_fn(shared$in_deg,  shared$in_connect)
  r_in_wt    <- rho_fn(shared$in_deg,  shared$in_degree)
  r_out_unwt <- rho_fn(shared$out_deg, shared$out_connect)
  r_out_wt   <- rho_fn(shared$out_deg, shared$out_degree)
  r_sevwod   <- rho_fn(shared$out_deg, shared$sevwod_combined)  # sevwod is out-degree based

  tibble(vp           = vp,
         rho_in_unwt  = r_in_unwt$rho,  p_in_unwt  = r_in_unwt$p,  rho_in_wt  = r_in_wt$rho,
         rho_out_unwt = r_out_unwt$rho, p_out_unwt = r_out_unwt$p, rho_out_wt = r_out_wt$rho,
         rho_sevwod   = r_sevwod$rho,
         n_shared     = nrow(shared))
})

# Print centrality summary with one-sample Wilcoxon per rho column;
# returns a summary tibble invisibly.
print_centrality_pilot <- function(res, label) {
  rho_cen_cols_local <- c("rho_in_unwt","rho_in_wt","rho_out_unwt","rho_out_wt","rho_sevwod")
  cat(sprintf("\n%s — Centrality Spearman (nar degree vs. PECAN):\n", label))
  summary_rows <- map_dfr(rho_cen_cols_local, function(col) {
    valid_rho <- res[[col]][!is.na(res[[col]])]
    cat(sprintf("  %-40s Mean rho: %.3f | SD: %.3f | Range: %s | n = %d\n", col,
                mean(valid_rho), sd(valid_rho), rng(res[[col]]), length(valid_rho)))
    row <- tibble(
      label    = label,
      rho_col  = col,
      n_valid  = length(valid_rho),
      mean_rho = round(mean(valid_rho), 3),
      sd_rho   = round(sd(valid_rho),   3),
      W        = NA_real_,
      p_value  = NA_real_
    )
    if (length(valid_rho) >= 3) {
      wt <- wilcox.test(valid_rho, mu = 0, alternative = "greater", exact = FALSE)
      cat(sprintf("  %-40s One-sample Wilcoxon (H1: rho > 0): V = %.0f, p = %.4f\n",
                  "", wt$statistic, wt$p.value))
      row$W       <- as.numeric(wt$statistic)
      row$p_value <- round(wt$p.value, 4)
    }
    row
  })
  invisible(summary_rows)
}

cat("\nCentrality Spearman (nar vs. PECAN, in- and out-degree separately):\n")
print(centrality_cor_results_pilot)
centrality_wilcox_pilot_unreg <- print_centrality_pilot(centrality_cor_results_pilot, "nar (unregulated)")

rho_cen_cols <- c("rho_in_unwt", "rho_in_wt", "rho_out_unwt", "rho_out_wt", "rho_sevwod")
write.csv(
  append_summary(centrality_cor_results_pilot, rho_cen_cols),
  "output/pilotstudy/tables/3c_centrality_spearman_unreg.csv", row.names = FALSE
)

# --- 3c (REGULATED). Node centrality comparison on regulated PECAN networks ---

nar_degree_reg_list_pilot <- list()

centrality_cor_results_reg_pilot <- map_dfr(vps_pilot, function(vp) {
  pn <- pecan_nodes_list_pilot[[vp]]
  pe <- pecan_edges_reg_list_pilot[[vp]]
  ne <- nar_edges_list_pilot[[vp]]

  na_row <- tibble(vp = vp,
                   rho_in_unwt  = NA_real_, p_in_unwt  = NA_real_, rho_in_wt  = NA_real_,
                   rho_out_unwt = NA_real_, p_out_unwt = NA_real_, rho_out_wt = NA_real_,
                   rho_sevwod   = NA_real_, n_shared   = NA_integer_)

  pecan_cen <- tryCatch(
    pecanCen(pe, pn, centrality_by = "width", sev_weighted = "value"),
    error = function(e) { warning(paste("[REG] pecanCen failed for", vp, ":", e$message)); NULL }
  )
  if (is.null(pecan_cen)) return(na_row)

  pecan_cen <- pecan_cen |> mutate(id = as.character(id))

  nar_out    <- ne |> count(from, name = "out_deg") |> rename(id = from)
  nar_in     <- ne |> count(to,   name = "in_deg")  |> rename(id = to)
  nar_degree <- full_join(nar_out, nar_in, by = "id") |>
    mutate(across(c(out_deg, in_deg), ~ replace_na(.x, 0L)),
           id = as.character(id))

  nar_degree_reg_list_pilot[[vp]] <<- nar_degree

  shared <- inner_join(
    nar_degree |> dplyr::select(id, in_deg, out_deg),
    pecan_cen  |> dplyr::select(id, in_connect, out_connect, in_degree, out_degree, sevwod_combined),
    by = "id"
  )

  if (nrow(shared) < 3) return(na_row |> mutate(n_shared = nrow(shared)))

  r_in_unwt  <- rho_fn(shared$in_deg,  shared$in_connect)
  r_in_wt    <- rho_fn(shared$in_deg,  shared$in_degree)
  r_out_unwt <- rho_fn(shared$out_deg, shared$out_connect)
  r_out_wt   <- rho_fn(shared$out_deg, shared$out_degree)
  r_sevwod   <- rho_fn(shared$out_deg, shared$sevwod_combined)

  tibble(vp           = vp,
         rho_in_unwt  = r_in_unwt$rho,  p_in_unwt  = r_in_unwt$p,  rho_in_wt  = r_in_wt$rho,
         rho_out_unwt = r_out_unwt$rho, p_out_unwt = r_out_unwt$p, rho_out_wt = r_out_wt$rho,
         rho_sevwod   = r_sevwod$rho,
         n_shared     = nrow(shared))
})

cat("\n[REGULATED] Centrality Spearman (nar vs. PECAN, in- and out-degree separately):\n")
print(centrality_cor_results_reg_pilot)
centrality_wilcox_pilot_reg <- print_centrality_pilot(centrality_cor_results_reg_pilot, "nar (regulated)")

write.csv(
  append_summary(centrality_cor_results_reg_pilot, rho_cen_cols),
  "output/pilotstudy/tables/3c_centrality_spearman_reg.csv", row.names = FALSE
)

# Combined unregulated + regulated centrality
centrality_cor_combined_pilot <- bind_rows(
  centrality_cor_results_pilot     |> mutate(pecan_version = "unregulated"),
  centrality_cor_results_reg_pilot |> mutate(pecan_version = "regulated")
) |> dplyr::select(pecan_version, everything())
write.csv(
  append_summary(centrality_cor_combined_pilot, rho_cen_cols, group_col = "pecan_version"),
  "output/pilotstudy/tables/3c_centrality_spearman_combined.csv", row.names = FALSE
)

centrality_wilcox_summary_pilot <- bind_rows(
  centrality_wilcox_pilot_unreg,
  centrality_wilcox_pilot_reg
)
write.csv(centrality_wilcox_summary_pilot,
          "output/pilotstudy/tables/3c_centrality_wilcoxon_summary.csv", row.names = FALSE)

# Shared-node degree vectors per VP — for tracing zero-variance exclusions.
# One row per VP × shared node (PECAN node set is the same for unreg and reg).
centrality_shared_nodes_pilot <- map_dfr(vps_pilot, function(vp) {
  shared_ids <- intersect(as.character(nar_degree_list_pilot[[vp]]$id),
                          as.character(pecan_nodes_list_pilot[[vp]]$id))
  nar_degree_list_pilot[[vp]] |>
    mutate(id = as.character(id)) |>
    filter(id %in% shared_ids) |>
    mutate(vp = vp)
}) |>
  dplyr::select(vp, id, in_deg, out_deg)

write.csv(centrality_shared_nodes_pilot,
          "output/pilotstudy/tables/3c_centrality_shared_node_vectors.csv",
          row.names = FALSE)

# --- 3c-TOP. Top-node agreement: does nar identify the same most-central
#             node as PECAN? (indegree and outdegree separately)
# For each participant we find the node(s) with the highest unweighted
# in-degree and out-degree in nar and in PECAN (over each network's full node
# set, not only shared nodes). Ties are handled conservatively: if multiple
# nodes share the maximum, we count agreement when at least one of the tied
# nar-top nodes is also a tied PECAN-top node.
#
# Three outcomes per VP per measure are distinguished:
#   match      — top nodes overlap across the two methods
#   no match   — top nodes differ, but the PECAN top node is present in nar
#   not shared — the top PECAN node does not appear in the nar network at all
#                (agreement is structurally impossible)
# ===========================================================================

# Helper: compute unweighted in/out degree for all nodes in one network
top_node_fn <- function(edges, nodes) {
  node_ids    <- as.character(nodes$id)
  node_labels <- setNames(as.character(nodes$label), node_ids)
  out_deg <- edges |> mutate(from = as.character(from)) |>
    filter(from %in% node_ids) |> count(from, name = "out_deg") |> rename(id = from)
  in_deg  <- edges |> mutate(to   = as.character(to))   |>
    filter(to   %in% node_ids) |> count(to,   name = "in_deg")  |> rename(id = to)
  tibble(id = node_ids) |>
    left_join(out_deg, by = "id") |>
    left_join(in_deg,  by = "id") |>
    mutate(out_deg = replace_na(out_deg, 0L),
           in_deg  = replace_na(in_deg,  0L),
           label   = node_labels[id])
}

# Helper: extract label string for a set of node IDs
get_lbl <- function(ids, nodes_df) {
  lbl <- setNames(as.character(nodes_df$label), as.character(nodes_df$id))
  paste(lbl[ids], collapse = " / ")
}

# Main function: per-VP top-node agreement (nar vs. one PECAN version)
# pecan_nodes_lst defaults to pecan_nodes_list_pilot (regulated PECAN shares same nodes)
top_node_agree_fn <- function(pecan_edges_lst, pecan_nodes_lst = pecan_nodes_list_pilot) {
  map_dfr(vps_pilot, function(vp) {
    nar_deg   <- top_node_fn(nar_edges_list_pilot[[vp]],  nar_nodes_list_pilot[[vp]])
    pecan_deg <- top_node_fn(pecan_edges_lst[[vp]],  pecan_nodes_lst[[vp]])

    # All nodes with the maximum degree (ties retained)
    top_nar_in    <- nar_deg   |> filter(in_deg  == max(in_deg))  |> pull(id)
    top_nar_out   <- nar_deg   |> filter(out_deg == max(out_deg)) |> pull(id)
    top_pecan_in  <- pecan_deg |> filter(in_deg  == max(in_deg))  |> pull(id)
    top_pecan_out <- pecan_deg |> filter(out_deg == max(out_deg)) |> pull(id)

    nar_ids <- as.character(nar_nodes_list_pilot[[vp]]$id)

    tibble(
      vp = vp,
      # --- Indegree ---
      top_nar_in_label    = get_lbl(top_nar_in,   nar_nodes_list_pilot[[vp]]),
      top_nar_in_deg      = max(nar_deg$in_deg),
      top_pecan_in_label  = get_lbl(top_pecan_in,  pecan_nodes_lst[[vp]]),
      top_pecan_in_deg    = max(pecan_deg$in_deg),
      pecan_top_in_shared = any(top_pecan_in %in% nar_ids),
      in_agree            = length(intersect(top_nar_in, top_pecan_in)) > 0,
      # --- Outdegree ---
      top_nar_out_label   = get_lbl(top_nar_out,  nar_nodes_list_pilot[[vp]]),
      top_nar_out_deg     = max(nar_deg$out_deg),
      top_pecan_out_label = get_lbl(top_pecan_out, pecan_nodes_lst[[vp]]),
      top_pecan_out_deg   = max(pecan_deg$out_deg),
      pecan_top_out_shared = any(top_pecan_out %in% nar_ids),
      out_agree            = length(intersect(top_nar_out, top_pecan_out)) > 0
    )
  })
}

# Summary helper: counts and conditional agreement rate
summarise_top_agree <- function(df, label) {
  n <- nrow(df)
  cat(sprintf("\n--- Top-node agreement: %s (N = %d) ---\n", label, n))
  for (measure in c("in", "out")) {
    shared_col <- paste0("pecan_top_", measure, "_shared")
    agree_col  <- paste0(measure, "_agree")
    n_shared   <- sum(df[[shared_col]])
    n_agree    <- sum(df[[agree_col]])
    n_agree_given_shared <- sum(df[[agree_col]] & df[[shared_col]])
    cat(sprintf(
      "  %sdegree: agree = %d/%d (%.0f%%) | PECAN top node in nar: %d/%d | agree|shared: %d/%d (%.0f%%)\n",
      measure, n_agree, n, 100 * n_agree / n,
      n_shared, n,
      n_agree_given_shared, n_shared,
      100 * n_agree_given_shared / max(n_shared, 1)
    ))
  }
}

# Compute for nar vs. unregulated and regulated PECAN
top_agree_nar_unreg_pilot <- top_node_agree_fn(pecan_edges_list_pilot)
top_agree_nar_reg_pilot   <- top_node_agree_fn(pecan_edges_reg_list_pilot)

summarise_top_agree(top_agree_nar_unreg_pilot, "nar vs. PECAN unregulated")
summarise_top_agree(top_agree_nar_reg_pilot,   "nar vs. PECAN regulated")

cat("\nPer-participant detail (nar vs. PECAN unregulated):\n")
print(as.data.frame(top_agree_nar_unreg_pilot))

cat("\nPer-participant detail (nar vs. PECAN regulated):\n")
print(as.data.frame(top_agree_nar_reg_pilot))

# Save per-participant tables
write.csv(top_agree_nar_unreg_pilot,
          "output/pilotstudy/tables/3c_top_node_agreement_nar_unreg.csv", row.names = FALSE)
write.csv(top_agree_nar_reg_pilot,
          "output/pilotstudy/tables/3c_top_node_agreement_nar_reg.csv",   row.names = FALSE)

# --- Overview summary table across both conditions ---
make_top_agree_overview <- function(df_list, labels) {
  map2_dfr(df_list, labels, function(df, label) {
    n <- nrow(df)
    tibble(
      comparison                  = label,
      N                           = n,
      # In-degree
      n_in_shared                 = sum(df$pecan_top_in_shared),
      pct_in_shared               = round(100 * n_in_shared / n, 1),
      n_in_agree                  = sum(df$in_agree),
      pct_in_agree                = round(100 * n_in_agree / n, 1),
      n_in_agree_given_shared     = sum(df$in_agree & df$pecan_top_in_shared),
      pct_in_agree_given_shared   = round(100 * n_in_agree_given_shared / max(n_in_shared, 1), 1),
      # Out-degree
      n_out_shared                = sum(df$pecan_top_out_shared),
      pct_out_shared              = round(100 * n_out_shared / n, 1),
      n_out_agree                 = sum(df$out_agree),
      pct_out_agree               = round(100 * n_out_agree / n, 1),
      n_out_agree_given_shared    = sum(df$out_agree & df$pecan_top_out_shared),
      pct_out_agree_given_shared  = round(100 * n_out_agree_given_shared / max(n_out_shared, 1), 1)
    )
  })
}

top_agree_overview_pilot <- make_top_agree_overview(
  list(top_agree_nar_unreg_pilot, top_agree_nar_reg_pilot),
  c("nar vs. PECAN (unregulated)", "nar vs. PECAN (regulated)")
)

cat("\n--- Top-node agreement overview ---\n")
print(top_agree_overview_pilot)

write.csv(top_agree_overview_pilot,
          "output/pilotstudy/tables/3c_top_node_agreement_overview.csv",
          row.names = FALSE)

cat("\nSaved top-node agreement tables.\n")

# ===========================================================================
# SECTION 4: DESCRIPTIVE & AGGREGATED DATA
# ===========================================================================
#relevant for the master thesis were the following output tables:
#4a_anova_network_metrics.csv #RM-ANOVA: tests descriptive network metrics for significant differences within a metric between network types
#4a_posthoc_network_metrics.csv #Bonferroni corrected paired t-tests: tests descriptive network metrics for significant differences within a metric between network types

cat("\n=== SECTION 4: DESCRIPTIVE & AGGREGATED DATA ===\n")

# Note: all_nar_nodes_pilot/edges and all_pecan_nodes_pilot/edges were pooled in Section 0.

# --- 4a. Network metrics compared across 3 network types (RM-ANOVA + post-hoc) ---
# Metrics: n nodes, n edges, n feedback loops (≤4), density
# Networks: nar, PECAN (unreg), PECAN (reg)
# Design: one-way within-subjects (repeated measures) ANOVA, 3 levels, n = 22 VPs.
# Sphericity: Greenhouse-Geisser correction applied automatically by afex.
# Post-hoc: all pairwise Bonferroni-corrected t-tests via emmeans::pairs().

metric_data_long_pilot <- map_dfr(vps_pilot, function(vp) {
  ne  <- distinct(nar_edges_list_pilot[[vp]], from, to)
  pe  <- pecan_edges_list_pilot[[vp]]
  pe_r <- pecan_edges_reg_list_pilot[[vp]]

  nn  <- nrow(nar_nodes_list_pilot[[vp]])
  np  <- nrow(pecan_nodes_list_pilot[[vp]])
  np_r <- nrow(pecan_nodes_list_pilot[[vp]])  # node set unchanged by regulation

  dens <- function(n_e, n_n) if (n_n > 1) n_e / (n_n * (n_n - 1)) else NA_real_

  loops_nar      <- count_cycles_le4(build_adj(ne,   nar_nodes_list_pilot[[vp]]$id))$total_le4
  loops_pecan    <- count_cycles_le4(build_adj(pe,   pecan_nodes_list_pilot[[vp]]$id))$total_le4
  loops_pecan_r  <- count_cycles_le4(build_adj(pe_r, pecan_nodes_list_pilot[[vp]]$id))$total_le4

  bind_rows(
    tibble(vp = vp, network = "nar",         n_nodes = nn,  n_edges = nrow(ne),   n_loops = loops_nar,     density = dens(nrow(ne),   nn)),
    tibble(vp = vp, network = "PECAN_unreg", n_nodes = np,  n_edges = nrow(pe),   n_loops = loops_pecan,   density = dens(nrow(pe),   np)),
    tibble(vp = vp, network = "PECAN_reg",   n_nodes = np_r, n_edges = nrow(pe_r), n_loops = loops_pecan_r, density = dens(nrow(pe_r), np_r))
  )
}) |>
  mutate(
    network = factor(network, levels = c("nar", "PECAN_unreg", "PECAN_reg")),
    vp      = factor(vp)
  )

write.csv(metric_data_long_pilot,
          "output/pilotstudy/tables/4a_network_metrics_long.csv", row.names = FALSE)

# Helper: RM-ANOVA + Bonferroni post-hoc for one metric
rm_anova_ph_fn <- function(metric_col, metric_label) {
  cat(sprintf("\n=== %s ===\n", metric_label))
  fit <- afex::aov_ez(
    id          = "vp",
    dv          = metric_col,
    data        = metric_data_long_pilot,
    within      = "network",
    anova_table = list(es = "pes", correction = "GG")
  )
  cat("RM-ANOVA (GG-corrected):\n"); print(fit)
  emm <- emmeans::emmeans(fit, ~ network)
  ph  <- pairs(emm, adjust = "bonferroni")
  cat("Post-hoc pairwise (Bonferroni):\n"); print(ph)
  list(anova = fit, posthoc = ph)
}

cat("\n--- 3c: RM-ANOVA — n nodes, n edges, n loops, density across network types ---\n")
anova_n_nodes_pilot <- rm_anova_ph_fn("n_nodes", "N Nodes")
anova_n_edges_pilot <- rm_anova_ph_fn("n_edges", "N Edges")
anova_n_loops_pilot <- rm_anova_ph_fn("n_loops", "N Feedback Loops (≤4)")
anova_density_pilot <- rm_anova_ph_fn("density", "Density")

anova_results_list_pilot <- list(
  n_nodes = anova_n_nodes_pilot, n_edges = anova_n_edges_pilot,
  n_loops = anova_n_loops_pilot, density = anova_density_pilot
)
anova_tbl_pilot <- map_dfr(names(anova_results_list_pilot), \(m)
  as.data.frame(anova_results_list_pilot[[m]]$anova$anova_table) |>
    mutate(metric = m, .before = 1)
)
posthoc_tbl_pilot <- map_dfr(names(anova_results_list_pilot), \(m)
  as.data.frame(summary(anova_results_list_pilot[[m]]$posthoc)) |>
    mutate(metric = m, .before = 1)
)
write.csv(anova_tbl_pilot,   "output/pilotstudy/tables/4a_anova_network_metrics.csv",   row.names = FALSE)
write.csv(posthoc_tbl_pilot, "output/pilotstudy/tables/4a_posthoc_network_metrics.csv", row.names = FALSE)
cat("[3c complete]\n")



# --- 4b. SIF: Symptom Importance Factor ---
# Formula: SIF = proportional_frequency × (mean_severity / max_scale)
# PECAN max scale = 10 | nar has no severity, so SIF_nar = proportional frequency only

sif_pecan_pilot <- all_pecan_nodes_pilot |>
  group_by(id, label) |>
  summarise(
    n_reported    = n(),
    mean_severity = round(mean(value, na.rm = TRUE), 2),
    sd_severity   = round(sd(value,   na.rm = TRUE), 2),
    min_severity  = min(value, na.rm = TRUE),
    max_severity  = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    prop_freq = round(n_reported / n_vps_pilot, 3),
    SIF       = round(prop_freq * (mean_severity / 10), 3)
  ) |>
  arrange(desc(SIF))

sif_nar_pilot <- all_nar_nodes_pilot |>
  group_by(id, label) |>
  summarise(n_reported = n(), .groups = "drop") |>
  mutate(
    prop_freq = round(n_reported / n_vps_pilot, 3),
    SIF_nar   = prop_freq  # no severity component available
  ) |>
  arrange(desc(SIF_nar))

cat("Top 10 symptoms by SIF (PECAN):\n")
print(head(sif_pecan_pilot, 10))
cat("\nTop 10 symptoms by SIF_nar:\n")
print(head(sif_nar_pilot, 10))

write.csv(sif_pecan_pilot, "output/pilotstudy/tables/4b_SIF_PECAN.csv", row.names = FALSE)
write.csv(sif_nar_pilot,   "output/pilotstudy/tables/4b_SIF_nar.csv",   row.names = FALSE)

# ---4c. EIF: Edge Impact Factor ---
# Formula: EIF = (n_VPs_with_edge / n_VPs_both_nodes_present) × (mean_width / max_scale)
# PECAN max_scale = 5 | nar all widths = 1 so EIF_nar = proportional frequency

# Pre-compute co-occurrence: for each VP, which node pairs co-occur?
pecan_cooccur_pilot <- build_cooccurrence(pecan_nodes_list_pilot, vps_pilot)
nar_cooccur_pilot   <- build_cooccurrence(nar_nodes_list_pilot,   vps_pilot)

pecan_cooccur_counts_pilot <- pecan_cooccur_pilot |>
  count(from, to, name = "n_both_present")
nar_cooccur_counts_pilot <- nar_cooccur_pilot |>
  count(from, to, name = "n_both_present")

# PECAN edge frequencies
pecan_edge_freq_pilot <- all_pecan_edges_pilot |>
  group_by(from, to) |>
  summarise(
    n_edge_reported = n(),
    mean_width      = round(mean(width, na.rm = TRUE), 2),
    sd_width        = round(sd(width,   na.rm = TRUE), 2),
    min_width       = min(width, na.rm = TRUE),
    max_width       = max(width, na.rm = TRUE),
    .groups = "drop"
  )

eif_pecan_pilot <- pecan_edge_freq_pilot |>
  left_join(pecan_cooccur_counts_pilot, by = c("from", "to")) |>
  mutate(
    prop_freq = if_else(
      !is.na(n_both_present) & n_both_present > 0,
      round(n_edge_reported / n_both_present, 3),
      NA_real_
    ),
    EIF = round(prop_freq * (mean_width / 5), 3)
  ) |>
  left_join(all_pecan_nodes_pilot |> dplyr::select(id, from_label = label) |> distinct(),
            by = c("from" = "id")) |>
  left_join(all_pecan_nodes_pilot |> dplyr::select(id, to_label = label) |> distinct(),
            by = c("to" = "id")) |>
  dplyr::select(from_label, to_label, from, to, n_edge_reported, n_both_present,
         prop_freq, mean_width, sd_width, min_width, max_width, EIF) |>
  arrange(desc(EIF))

# nar edge frequencies
nar_edge_freq_pilot <- all_nar_edges_pilot |>
  distinct(from, to, vp) |>  # deduplicate within VP
  group_by(from, to) |>
  summarise(n_edge_reported = n(), .groups = "drop")

eif_nar_pilot <- nar_edge_freq_pilot |>
  left_join(nar_cooccur_counts_pilot, by = c("from", "to")) |>
  mutate(
    prop_freq = if_else(
      !is.na(n_both_present) & n_both_present > 0,
      round(n_edge_reported / n_both_present, 3),
      NA_real_
    ),
    EIF_nar = prop_freq
  ) |>
  left_join(all_nar_nodes_pilot |> dplyr::select(id, from_label = label) |> distinct(),
            by = c("from" = "id")) |>
  left_join(all_nar_nodes_pilot |> dplyr::select(id, to_label = label) |> distinct(),
            by = c("to" = "id")) |>
  dplyr::select(from_label, to_label, from, to, n_edge_reported, n_both_present,
         prop_freq, EIF_nar) |>
  arrange(desc(EIF_nar))

cat("Top 10 edges by EIF (PECAN):\n")
print(head(eif_pecan_pilot, 10))
cat("\nTop 10 edges by EIF_nar:\n")
print(head(eif_nar_pilot, 10))

write.csv(eif_pecan_pilot, "output/pilotstudy/tables/4c_EIF_PECAN_unreg.csv", row.names = FALSE)
write.csv(eif_nar_pilot,   "output/pilotstudy/tables/4c_EIF_nar.csv",         row.names = FALSE)

# --- 4c (REGULATED). EIF on regulated PECAN networks ---
# Co-occurrence counts are unchanged (node-based); only edge frequencies differ.

pecan_edge_freq_reg_pilot <- all_pecan_edges_reg_pilot |>
  group_by(from, to) |>
  summarise(
    n_edge_reported = n(),
    mean_width      = round(mean(width, na.rm = TRUE), 2),
    sd_width        = round(sd(width,   na.rm = TRUE), 2),
    min_width       = min(width, na.rm = TRUE),
    max_width       = max(width, na.rm = TRUE),
    .groups = "drop"
  )

eif_pecan_reg_pilot <- pecan_edge_freq_reg_pilot |>
  left_join(pecan_cooccur_counts_pilot, by = c("from", "to")) |>
  mutate(
    prop_freq = if_else(
      !is.na(n_both_present) & n_both_present > 0,
      round(n_edge_reported / n_both_present, 3),
      NA_real_
    ),
    EIF = round(prop_freq * (mean_width / 5), 3)
  ) |>
  left_join(all_pecan_nodes_pilot |> dplyr::select(id, from_label = label) |> distinct(),
            by = c("from" = "id")) |>
  left_join(all_pecan_nodes_pilot |> dplyr::select(id, to_label = label) |> distinct(),
            by = c("to" = "id")) |>
  dplyr::select(from_label, to_label, from, to, n_edge_reported, n_both_present,
         prop_freq, mean_width, sd_width, min_width, max_width, EIF) |>
  arrange(desc(EIF))

cat("\n[REGULATED] Top 10 edges by EIF (PECAN):\n")
print(head(eif_pecan_reg_pilot, 10))

write.csv(eif_pecan_reg_pilot, "output/pilotstudy/tables/4c_EIF_PECAN_reg.csv", row.names = FALSE)

# --- 4d. Feedback Loops (matrix algebra, cycles up to length 4) ---
# Count directed simple cycles of length 2, 3, and 4 using trace of adjacency
# matrix powers. This is O(n³) and bounded regardless of network size.
#
# For binary directed adjacency matrix A (no self-loops):
#   n_2cycles = trace(A²) / 2
#   n_3cycles = trace(A³) / 3
#   n_4cycles = (trace(A⁴) − 2·Σᵢ(A²[i,i])² + trace(A²)) / 4
#
# The n_4cycles correction subtracts non-simple length-4 closed walks:
#   Case v₁=v₃ (wings node revisited): contributes Σᵢ (A²[i,i])² to trace(A⁴)
#   Case v₂=v₀ (start revisited mid-walk): contributes Σᵢ (A²[i,i])² to trace(A⁴)
#   Overlap (both cases simultaneously, e.g. i→j→i→j→i): contributes trace(A²)
# Inclusion-exclusion: 2·Σ(A²diag²) − trace(A²), subtracted from trace(A⁴).

# PECAN loop counts
loop_pecan_pilot <- map_dfr(vps_pilot, function(vp) {
  adj <- build_adj(pecan_edges_list_pilot[[vp]], pecan_nodes_list_pilot[[vp]]$id)
  cyc <- count_cycles_le4(adj)
  tibble(vp = vp,
         pecan_n2  = cyc$n2,  pecan_n3  = cyc$n3,
         pecan_n4  = cyc$n4,  pecan_le4 = cyc$total_le4)
})

# nar loop counts (edges are binary; deduplicate within VP)
loop_nar_pilot <- map_dfr(vps_pilot, function(vp) {
  adj <- build_adj(nar_edges_list_pilot[[vp]], nar_nodes_list_pilot[[vp]]$id)
  cyc <- count_cycles_le4(adj)
  tibble(vp = vp,
         nar_n2  = cyc$n2,  nar_n3  = cyc$n3,
         nar_n4  = cyc$n4,  nar_le4 = cyc$total_le4)
})

# Regulated PECAN loop counts
loop_pecan_reg_pilot <- map_dfr(vps_pilot, function(vp) {
  adj <- build_adj(pecan_edges_reg_list_pilot[[vp]], pecan_nodes_list_pilot[[vp]]$id)
  cyc <- count_cycles_le4(adj)
  tibble(vp = vp,
         pecan_reg_n2  = cyc$n2,  pecan_reg_n3  = cyc$n3,
         pecan_reg_n4  = cyc$n4,  pecan_reg_le4 = cyc$total_le4)
})

loop_results_pilot <- left_join(loop_pecan_pilot, loop_nar_pilot, by = "vp") |>
  left_join(loop_pecan_reg_pilot, by = "vp")

cat("\nFeedback loops per participant (unregulated PECAN | nar | regulated PECAN, cycles ≤ 4):\n")
print(loop_results_pilot)

cat("\n--- Summary ---\n")
cat("PECAN unregulated  mean loops (≤4):", round(mean(loop_results_pilot$pecan_le4),     1),
    "| SD:", round(sd(loop_results_pilot$pecan_le4),     1), "\n")
cat("PECAN regulated    mean loops (≤4):", round(mean(loop_results_pilot$pecan_reg_le4), 1),
    "| SD:", round(sd(loop_results_pilot$pecan_reg_le4), 1), "\n")
cat("nar                mean loops (≤4):", round(mean(loop_results_pilot$nar_le4),       1),
    "| SD:", round(sd(loop_results_pilot$nar_le4),       1), "\n")

# Paired test: nar vs. PECAN unregulated
loop_diff_unreg_pilot <- loop_results_pilot$nar_le4 - loop_results_pilot$pecan_le4
sw_loops_unreg_pilot  <- shapiro.test(loop_diff_unreg_pilot)
cat("\nShapiro-Wilk (nar − PECAN unreg diff): W =", round(sw_loops_unreg_pilot$statistic, 3),
    ", p =", round(sw_loops_unreg_pilot$p.value, 4), "\n")
cat("\nPaired test (nar vs. PECAN unregulated, n_loops_le4):\n")
if (sw_loops_unreg_pilot$p.value > 0.05) {
  print(t.test(loop_results_pilot$nar_le4, loop_results_pilot$pecan_le4, paired = TRUE))
} else {
  print(wilcox.test(loop_results_pilot$nar_le4, loop_results_pilot$pecan_le4, paired = TRUE))
}

# Paired test: nar vs. PECAN regulated
loop_diff_reg_pilot  <- loop_results_pilot$nar_le4 - loop_results_pilot$pecan_reg_le4
sw_loops_reg_pilot   <- shapiro.test(loop_diff_reg_pilot)
cat("\nShapiro-Wilk (nar − PECAN reg diff): W =", round(sw_loops_reg_pilot$statistic, 3),
    ", p =", round(sw_loops_reg_pilot$p.value, 4), "\n")
cat("\n[REGULATED] Paired test (nar vs. PECAN regulated, n_loops_le4):\n")
if (sw_loops_reg_pilot$p.value > 0.05) {
  print(t.test(loop_results_pilot$nar_le4, loop_results_pilot$pecan_reg_le4, paired = TRUE))
} else {
  print(wilcox.test(loop_results_pilot$nar_le4, loop_results_pilot$pecan_reg_le4, paired = TRUE))
}

write.csv(loop_results_pilot, "output/pilotstudy/tables/4d_feedback_loops.csv", row.names = FALSE)

# --- 4d-ii. Common loops across participants ---
# Enumerate all simple directed cycles up to length 4 by node label,
# then count how many participants share each cycle.
#
# Canonical form: lexicographically smallest cyclic rotation of the label
# sequence (direction preserved; only the starting node is shifted).



# --- PECAN common loops ---
pecan_cycle_list_pilot <- setNames(
  lapply(vps_pilot, function(vp)
    enumerate_cycles_labeled(pecan_edges_list_pilot[[vp]], pecan_nodes_list_pilot[[vp]])),
  vps_pilot
)

pecan_common_loops_pilot <- pecan_cycle_list_pilot |>
  imap_dfr(~ tibble(vp = .y, cycle = .x)) |>
  filter(!is.na(cycle)) |>
  group_by(cycle) |>
  summarise(
    n_participants = n(),
    participants   = paste(sort(vp), collapse = ", "),
    .groups = "drop"
  ) |>
  mutate(cycle_length = str_count(cycle, "\u2192") + 1L) |>
  arrange(desc(n_participants), cycle_length, cycle)

cat("\nCommon feedback loops across participants (PECAN, cycles ≤ 4):\n")
cat("Total unique cycles found:", nrow(pecan_common_loops_pilot), "\n")
cat("Cycles shared by ≥ 2 participants:\n")
print(pecan_common_loops_pilot |> filter(n_participants >= 2))

# --- nar common loops ---
nar_cycle_list_pilot <- setNames(
  lapply(vps_pilot, function(vp)
    enumerate_cycles_labeled(nar_edges_list_pilot[[vp]], nar_nodes_list_pilot[[vp]])),
  vps_pilot
)

nar_common_loops_pilot <- nar_cycle_list_pilot |>
  imap_dfr(~ tibble(vp = .y, cycle = .x)) |>
  filter(!is.na(cycle)) |>
  group_by(cycle) |>
  summarise(
    n_participants = n(),
    participants   = paste(sort(vp), collapse = ", "),
    .groups = "drop"
  ) |>
  mutate(cycle_length = str_count(cycle, "\u2192") + 1L) |>
  arrange(desc(n_participants), cycle_length, cycle)

cat("\nCommon feedback loops across participants (nar, cycles ≤ 4):\n")
cat("Total unique cycles found:", nrow(nar_common_loops_pilot), "\n")
cat("Cycles shared by ≥ 2 participants:\n")
print(nar_common_loops_pilot |> filter(n_participants >= 2))

write.csv(pecan_common_loops_pilot, "output/pilotstudy/tables/4d_common_loops_PECAN.csv", row.names = FALSE)
write.csv(nar_common_loops_pilot,   "output/pilotstudy/tables/4d_common_loops_nar.csv",   row.names = FALSE)

# --- 4e. Aggregated Group Networks ---
# Inclusion threshold_pilot: ≥ 10% of participants = ≥ 3 of 22

threshold_pilot <- ceiling(0.10 * n_vps_pilot)
cat("\nAggregation threshold_pilot:", threshold_pilot, "of", n_vps_pilot, "participants\n")

# PECAN group network
pecan_node_include_pilot <- sif_pecan_pilot |> filter(n_reported >= threshold_pilot)
pecan_edge_include_pilot <- eif_pecan_pilot |>
  filter(n_edge_reported >= threshold_pilot,
         from %in% pecan_node_include_pilot$id,
         to   %in% pecan_node_include_pilot$id)

cat("PECAN group network: ", nrow(pecan_node_include_pilot), "nodes,",
    nrow(pecan_edge_include_pilot), "edges\n")

if (nrow(pecan_edge_include_pilot) > 0 && nrow(pecan_node_include_pilot) > 0) {
  grp_ids  <- unique(c(pecan_edge_include_pilot$from, pecan_edge_include_pilot$to))
  grp_info_pilot <- pecan_node_include_pilot |>
    filter(id %in% grp_ids) |>
    arrange(match(id, grp_ids))

  pecan_grp_adj_pilot <- matrix(
    0, length(grp_ids), length(grp_ids),
    dimnames = list(as.character(grp_ids), as.character(grp_ids))
  )
  for (i in seq_len(nrow(pecan_edge_include_pilot))) {
    f  <- as.character(pecan_edge_include_pilot$from[i])
    t_ <- as.character(pecan_edge_include_pilot$to[i])
    if (f %in% rownames(pecan_grp_adj_pilot) && t_ %in% colnames(pecan_grp_adj_pilot)) {
      pecan_grp_adj_pilot[f, t_] <- pecan_edge_include_pilot$mean_width[i]
    }
  }

  png("output/pilotstudy/plots/group_network_PECAN_unreg.png", width = 1600, height = 1200, res = 150)
  qgraph(pecan_grp_adj_pilot,
         labels     = grp_info_pilot$label,
         edge.color = "steelblue",
         title      = paste0("PECAN Group Network UNREGULATED (≥", threshold_pilot, "/", n_vps_pilot, " VPs)"),
         directed   = TRUE,
         layout     = "spring")
  dev.off()
  cat("Saved: output/plots/group_network_PECAN_unreg.png\n")
}

# PECAN group network (REGULATED)
pecan_node_include_reg_pilot <- sif_pecan_pilot |> filter(n_reported >= threshold_pilot)
pecan_edge_include_reg_pilot <- eif_pecan_reg_pilot |>
  filter(n_edge_reported >= threshold_pilot,
         from %in% pecan_node_include_reg_pilot$id,
         to   %in% pecan_node_include_reg_pilot$id)

cat("[REGULATED] PECAN group network: ", nrow(pecan_node_include_reg_pilot), "nodes,",
    nrow(pecan_edge_include_reg_pilot), "edges\n")

if (nrow(pecan_edge_include_reg_pilot) > 0 && nrow(pecan_node_include_reg_pilot) > 0) {
  grp_ids_reg  <- unique(c(pecan_edge_include_reg_pilot$from, pecan_edge_include_reg_pilot$to))
  grp_info_reg_pilot <- pecan_node_include_reg_pilot |>
    filter(id %in% grp_ids_reg) |>
    arrange(match(id, grp_ids_reg))

  pecan_grp_adj_reg_pilot <- matrix(
    0, length(grp_ids_reg), length(grp_ids_reg),
    dimnames = list(as.character(grp_ids_reg), as.character(grp_ids_reg))
  )
  for (i in seq_len(nrow(pecan_edge_include_reg_pilot))) {
    f  <- as.character(pecan_edge_include_reg_pilot$from[i])
    t_ <- as.character(pecan_edge_include_reg_pilot$to[i])
    if (f %in% rownames(pecan_grp_adj_reg_pilot) && t_ %in% colnames(pecan_grp_adj_reg_pilot)) {
      pecan_grp_adj_reg_pilot[f, t_] <- pecan_edge_include_reg_pilot$mean_width[i]
    }
  }

  png("output/pilotstudy/plots/group_network_PECAN_reg.png", width = 1600, height = 1200, res = 150)
  qgraph(pecan_grp_adj_reg_pilot,
         labels     = grp_info_reg_pilot$label,
         edge.color = "steelblue",
         title      = paste0("PECAN Group Network REGULATED (≥", threshold_pilot, "/", n_vps_pilot, " VPs)"),
         directed   = TRUE,
         layout     = "spring")
  dev.off()
  cat("Saved: output/plots/group_network_PECAN_reg.png\n")
}

# nar group network
nar_node_include_pilot <- sif_nar_pilot |> filter(n_reported >= threshold_pilot)
nar_edge_include_pilot <- eif_nar_pilot |>
  filter(n_edge_reported >= threshold_pilot,
         from %in% nar_node_include_pilot$id,
         to   %in% nar_node_include_pilot$id)

cat("nar group network: ", nrow(nar_node_include_pilot), "nodes,",
    nrow(nar_edge_include_pilot), "edges\n")

if (nrow(nar_edge_include_pilot) > 0 && nrow(nar_node_include_pilot) > 0) {
  grp_ids_nar  <- unique(c(nar_edge_include_pilot$from, nar_edge_include_pilot$to))
  grp_info_nar_pilot <- nar_node_include_pilot |> filter(id %in% grp_ids_nar)

  nar_grp_adj_pilot <- matrix(
    0, length(grp_ids_nar), length(grp_ids_nar),
    dimnames = list(as.character(grp_ids_nar), as.character(grp_ids_nar))
  )
  for (i in seq_len(nrow(nar_edge_include_pilot))) {
    f  <- as.character(nar_edge_include_pilot$from[i])
    t_ <- as.character(nar_edge_include_pilot$to[i])
    if (f %in% rownames(nar_grp_adj_pilot) && t_ %in% colnames(nar_grp_adj_pilot)) {
      nar_grp_adj_pilot[f, t_] <- 1
    }
  }

  png("output/pilotstudy/plots/group_network_nar.png", width = 1600, height = 1200, res = 150)
  qgraph(nar_grp_adj_pilot,
         labels     = grp_info_nar_pilot$label,
         edge.color = "darkgreen",
         title      = paste0("nar Group Network (≥", threshold_pilot, "/", n_vps_pilot, " VPs)"),
         directed   = TRUE,
         layout     = "spring")
  dev.off()
  cat("Saved: output/plots/group_network_nar.png\n")
}

# --- 4f. Symptom frequency table (combined nar + PECAN) ---

symptom_table_pilot <- full_join(
  sif_nar_pilot   |> dplyr::select(id, label, n_nar = n_reported, pct_nar = prop_freq),
  sif_pecan_pilot |> dplyr::select(id, pecan_label = label, n_pecan = n_reported,
                      pct_pecan = prop_freq, mean_severity, sd_severity,
                      min_severity, max_severity, SIF),
  by = "id"
) |>
  mutate(
    label = coalesce(label, pecan_label),
    n_nar   = replace_na(n_nar, 0L),
    n_pecan = replace_na(n_pecan, 0L)
  ) |>
  dplyr::select(-pecan_label) |>
  arrange(desc(n_pecan + n_nar))

cat("\nSymptom frequency table (top 20):\n")
print(head(symptom_table_pilot, 20))
write.csv(symptom_table_pilot, "output/pilotstudy/tables/4f_symptom_frequency_table.csv", row.names = FALSE)

# --- 4g. Edge frequency table (PECAN with descriptive statistics) ---

cat("\nTop 20 edges by frequency (PECAN unregulated):\n")
print(head(eif_pecan_pilot |> arrange(desc(n_edge_reported)), 20))
write.csv(
  eif_pecan_pilot |> arrange(desc(n_edge_reported)),
  "output/pilotstudy/tables/4g_edge_frequency_PECAN_unreg.csv", row.names = FALSE
)

cat("\n[REGULATED] Top 20 edges by frequency (PECAN regulated):\n")
print(head(eif_pecan_reg_pilot |> arrange(desc(n_edge_reported)), 20))
write.csv(
  eif_pecan_reg_pilot |> arrange(desc(n_edge_reported)),
  "output/pilotstudy/tables/4g_edge_frequency_PECAN_reg.csv", row.names = FALSE
)

# --- 4h. nar vs. PECAN edge comparison ---
# How often was each directed edge reported in nar vs. PECAN?

# Label lookup combining both networks
label_lookup_pilot <- bind_rows(
  all_nar_nodes_pilot   |> dplyr::select(id, label),
  all_pecan_nodes_pilot |> dplyr::select(id, label)
) |>
  distinct(id, .keep_all = TRUE)

edge_comparison_pilot <- full_join(
  nar_edge_freq_pilot   |> dplyr::select(from, to, n_nar = n_edge_reported),
  pecan_edge_freq_pilot |> dplyr::select(from, to, n_pecan = n_edge_reported),  by = c("from", "to")
) |>
  mutate(
    n_nar   = replace_na(n_nar, 0L),
    n_pecan = replace_na(n_pecan, 0L),
    diff    = n_nar - n_pecan,
    more_in = case_when(
      diff > 0 ~ "nar",
      diff < 0 ~ "PECAN",
      TRUE     ~ "equal"
    )
  ) |>
  left_join(label_lookup_pilot |> rename(from_label = label), by = c("from" = "id")) |>
  left_join(label_lookup_pilot |> rename(to_label   = label), by = c("to"   = "id")) |>
  dplyr::select(from_label, to_label, from, to, n_nar, n_pecan, diff, more_in) |>
  arrange(desc(abs(diff)))

cat("\nEdge comparison — most divergent (nar vs. PECAN):\n")
print(head(edge_comparison_pilot, 20))
cat("\nSummary of edge direction:\n")
print(count(edge_comparison_pilot, more_in))

write.csv(edge_comparison_pilot, "output/pilotstudy/tables/4h_edge_comparison_nar_vs_PECAN_unreg.csv",
          row.names = FALSE)

# --- 4h (REGULATED). nar vs. PECAN regulated edge comparison ---

edge_comparison_reg_pilot <- full_join(
  nar_edge_freq_pilot        |> dplyr::select(from, to, n_nar   = n_edge_reported),
  pecan_edge_freq_reg_pilot  |> dplyr::select(from, to, n_pecan = n_edge_reported),
  by = c("from", "to")
) |>
  mutate(
    n_nar   = replace_na(n_nar, 0L),
    n_pecan = replace_na(n_pecan, 0L),
    diff    = n_nar - n_pecan,
    more_in = case_when(
      diff > 0 ~ "nar",
      diff < 0 ~ "PECAN_reg",
      TRUE     ~ "equal"
    )
  ) |>
  left_join(label_lookup_pilot |> rename(from_label = label), by = c("from" = "id")) |>
  left_join(label_lookup_pilot |> rename(to_label   = label), by = c("to"   = "id")) |>
  dplyr::select(from_label, to_label, from, to, n_nar, n_pecan, diff, more_in) |>
  arrange(desc(abs(diff)))

cat("\n[REGULATED] Edge comparison — most divergent (nar vs. PECAN regulated):\n")
print(head(edge_comparison_reg_pilot, 20))
cat("\n[REGULATED] Summary of edge direction:\n")
print(count(edge_comparison_reg_pilot, more_in))

write.csv(edge_comparison_reg_pilot, "output/pilotstudy/tables/4h_edge_comparison_nar_vs_PECAN_reg.csv",
          row.names = FALSE)

# ===========================================================================
# SECTION 5: SUBJECTIVE REPRESENTATIVENESS (RM-ANOVA)
# ===========================================================================
#relevant for the master thesis were the following output tables and plots:
#5_cronbach_alpha_rep_nar.csv #internal consistency of the representativeness and feasibility scale
#5_cronbach_alpha_rep_PECAN.csv #internal consistency of the representativeness and feasibility scale
#5b_representativeness_descriptives_reduced.csv #descriptives of subjective representativeness and feasibility ratings
#5b_representativeness_anova_reduced.csv #ANOVA comparing the subjective representativeness of the networks (nar_full, nar_ex, PECAN reg)

# Participants rated representativeness of their networks after viewing them.
# NW02 = PECAN, NW07 = nar_ex (NW06/nar_im not collected in pilot study)
#
# Reversed items (higher score = LESS representative → flip before summing):
#   NW02_14, NW02_15, NW07_14, NW07_15
# Reversal formula: 6 - x  (scale 1–5)

# --- 5.1 T2 data ---
# t2_pilot was loaded in Section 2a-T2; no reload needed here.
cat("\nT2 data ready:", nrow(t2_pilot), "participants\n")

# --- 5.2 Item lists and reversal ---

nw02_items_pilot <- names(t2_pilot)[grepl("^NW02_", names(t2_pilot))]
nw07_items_pilot <- names(t2_pilot)[grepl("^NW07_", names(t2_pilot))]
rev_items  <- c("NW02_14", "NW02_15", "NW07_14", "NW07_15")

cat("NW02 items (PECAN, n =", length(nw02_items_pilot), "):", paste(nw02_items_pilot, collapse = ", "), "\n")
cat("NW07 items (nar_ex, n =", length(nw07_items_pilot), "):", paste(nw07_items_pilot, collapse = ", "), "\n")
cat("Reversed:   ", paste(rev_items, collapse = ", "), "\n")

t2_scored_pilot <- t2_pilot |>
  mutate(across(all_of(rev_items), ~ 6 - .x))

# --- 5.3 Cronbach's alpha for representativeness questionnaire ---

alpha_rep_pecan_pilot <- withCallingHandlers(
  psych::alpha(t2_scored_pilot |> dplyr::select(all_of(nw02_items_pilot)), check.keys = FALSE),
  warning = function(w) invokeRestart("muffleWarning")
)
alpha_rep_nar_pilot <- withCallingHandlers(
  psych::alpha(t2_scored_pilot |> dplyr::select(all_of(nw07_items_pilot)), check.keys = FALSE),
  warning = function(w) invokeRestart("muffleWarning")
)

cat("\nCronbach's Alpha — Representativeness questionnaire (PECAN / NW02):\n")
print(alpha_rep_pecan_pilot$total)
cat("\nCronbach's Alpha — Representativeness questionnaire (nar_ex / NW07):\n")
print(alpha_rep_nar_pilot$total)

write.csv(as.data.frame(alpha_rep_pecan_pilot$total),
          "output/pilotstudy/tables/5_cronbach_alpha_rep_PECAN.csv", row.names = FALSE)
write.csv(as.data.frame(alpha_rep_nar_pilot$total),
          "output/pilotstudy/tables/5_cronbach_alpha_rep_nar.csv",   row.names = FALSE)

# --- 5.4 Sum scores ---

rep_data_pilot <- t2_scored_pilot |>
  dplyr::select(VP, all_of(c(nw02_items_pilot, nw07_items_pilot))) |>
  mutate(
    PECAN  = rowSums(across(all_of(nw02_items_pilot)), na.rm = FALSE),
    nar_ex = rowSums(across(all_of(nw07_items_pilot)), na.rm = FALSE)
  ) |>
  dplyr::select(VP, PECAN, nar_ex)

cat("\nSum scores per participant:\n")
print(rep_data_pilot)

cat("\nMissing values: PECAN =", sum(is.na(rep_data_pilot$PECAN)),
    "| nar_ex =", sum(is.na(rep_data_pilot$nar_ex)), "\n")

rep_complete_pilot <- rep_data_pilot |> filter(!is.na(PECAN) & !is.na(nar_ex))
cat("Complete cases for ANOVA: n =", nrow(rep_complete_pilot), "\n")

# --- 5.5 Descriptives ---

rep_long_pilot <- rep_complete_pilot |>
  pivot_longer(cols = c(PECAN, nar_ex),
               names_to  = "model",
               values_to = "representativeness") |>
  mutate(VP    = factor(VP),
         model = factor(model, levels = c("PECAN", "nar_ex")))

max_score_pilot <- length(nw02_items_pilot) * 5

rep_desc_pilot <- rep_long_pilot |>
  group_by(model) |>
  summarise(M   = round(mean(representativeness), 2),
            SD  = round(sd(representativeness),   2),
            Min = min(representativeness),
            Max = max(representativeness),
            .groups = "drop") |>
  mutate(
    M_pct   = round(M   / max_score_pilot * 100, 1),
    SD_pct  = round(SD  / max_score_pilot * 100, 1),
    Min_pct = round(Min / max_score_pilot * 100, 1),
    Max_pct = round(Max / max_score_pilot * 100, 1)
  )

cat("\nDescriptives (sum scores, max =", max_score_pilot, "| percentages of max):\n")
print(rep_desc_pilot)

# --- 5.6 Assumption check: normality of difference scores ---

diff_scores_pilot <- rep_complete_pilot$PECAN - rep_complete_pilot$nar_ex
sw_rep_pilot      <- shapiro.test(diff_scores_pilot)
cat("\nShapiro-Wilk (PECAN - nar_ex difference): W =",
    round(sw_rep_pilot$statistic, 3), ", p =", round(sw_rep_pilot$p.value, 4), "\n")

# --- 5.7 RM ANOVA ---
# Note: with 2 conditions, sphericity is trivially met (no Mauchly's test needed)
# and RM ANOVA is mathematically equivalent to a paired t-test.

rep_anova_pilot <- aov_ez(
  id      = "VP",
  dv      = "representativeness",
  data    = rep_long_pilot,
  within  = "model",
  anova_table = list(es = "pes")
)

cat("\nRepeated Measures ANOVA — subjective representativeness (PECAN vs. nar_ex):\n")
print(rep_anova_pilot)

# --- 5.8 Paired t-test (equivalent with 2 conditions; easier to report) ---

cat("\nPaired t-test (PECAN vs. nar_ex):\n")
print(t.test(rep_complete_pilot$PECAN, rep_complete_pilot$nar_ex, paired = TRUE))

# --- 5.8b Wilcoxon signed-rank (sensitivity check) ---
# Difference scores are non-normally distributed (Shapiro-Wilk p = .017),
# so the Wilcoxon signed-rank test is reported alongside the parametric results.

cat("\nWilcoxon signed-rank test (sensitivity check):\n")
print(wilcox.test(rep_complete_pilot$PECAN, rep_complete_pilot$nar_ex, paired = TRUE, exact = FALSE))

# --- 5.8c Paired plot ---

ggplot(rep_long_pilot, aes(x = model, y = representativeness)) +
  geom_line(aes(group = VP), colour = "grey60", linewidth = 0.5) +
  geom_point(aes(group = VP), size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, colour = "steelblue") +
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               colour = "steelblue", linewidth = 1.2, linetype = "dashed") +
  scale_x_discrete(labels = c("PECAN", "nar_ex")) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
  labs(
    x       = "Network model",
    y       = "Representativeness (sum score)",
    title   = "Subjective representativeness: PECAN vs. nar_ex",
    caption = paste0("Grey lines = individual participants (n = ", nrow(rep_complete_pilot),
                     "). Blue diamond = group mean.")
  )

ggsave("output/pilotstudy/plots/representativeness_paired_plot.png",
       width = 6, height = 5, dpi = 150)

# --- 5.9 Save outputs ---

write.csv(rep_data_pilot,
          "output/pilotstudy/tables/5_representativeness_sumscores.csv", row.names = FALSE)
write.csv(rep_desc_pilot,
          "output/pilotstudy/tables/5_representativeness_descriptives.csv", row.names = FALSE)
write.csv(as.data.frame(rep_anova_pilot$anova_table),
          "output/pilotstudy/tables/5_representativeness_anova.csv", row.names = TRUE)

cat("Saved representativeness tables.\n")

# ---------------------------------------------------------------------------
# SECTION 5b: SUBJECTIVE REPRESENTATIVENESS — NETWORK-SPECIFIC ITEMS ONLY
# ---------------------------------------------------------------------------
# Items NW02_05–NW02_09 (and the equivalent NW07_05–NW07_09) ask about
# networks in general rather than the specific network shown. They are
# excluded here to obtain a sum score reflecting only network-specific
# representativeness judgements.
#
# Everything else is identical to Section 5 (same reversal, same ANOVA).
# Variables carry the suffix _red ("reduced item set").
# Reversed items: NW*_14 and NW*_15 are already reversed in t2_scored_pilot.
# Note: NW06 (nar_im) was not collected in the pilot — only PECAN vs. nar_ex.

# --- 5b.2 Reduced item lists (exclude items 05–09) ---

general_items <- paste0("0", 5:9)   # "05","06","07","08","09"

nw02_items_red_pilot <- nw02_items_pilot[!sub("^NW02_", "", nw02_items_pilot) %in% general_items]
nw07_items_red_pilot <- nw07_items_pilot[!sub("^NW07_", "", nw07_items_pilot) %in% general_items]

cat("\n[5b] Network-specific items only (NW*_05–NW*_09 excluded):\n")
cat("NW02 items (PECAN,  n =", length(nw02_items_red_pilot), "):", paste(nw02_items_red_pilot, collapse = ", "), "\n")
cat("NW07 items (nar_ex, n =", length(nw07_items_red_pilot), "):", paste(nw07_items_red_pilot, collapse = ", "), "\n")
cat("Excluded (general network items):",
    paste(paste0("NW02_", general_items), collapse = ", "), "and NW07 equivalents\n")

# --- 5b.3 Sum scores (reduced) ---

rep_data_red_pilot <- t2_scored_pilot |>
  dplyr::select(VP, all_of(c(nw02_items_red_pilot, nw07_items_red_pilot))) |>
  mutate(
    PECAN  = rowSums(across(all_of(nw02_items_red_pilot)), na.rm = FALSE),
    nar_ex = rowSums(across(all_of(nw07_items_red_pilot)), na.rm = FALSE)
  ) |>
  dplyr::select(VP, PECAN, nar_ex)

cat("\n[5b] Reduced sum scores per participant:\n")
print(rep_data_red_pilot)
cat("\nMissing values: PECAN =", sum(is.na(rep_data_red_pilot$PECAN)),
    "| nar_ex =", sum(is.na(rep_data_red_pilot$nar_ex)), "\n")

rep_complete_red_pilot <- rep_data_red_pilot |> filter(!is.na(PECAN) & !is.na(nar_ex))
cat("Complete cases for RM-ANOVA: n =", nrow(rep_complete_red_pilot), "\n")

# --- 5b.4 Descriptives ---

rep_long_red_pilot <- rep_complete_red_pilot |>
  pivot_longer(cols = c(PECAN, nar_ex),
               names_to  = "model",
               values_to = "representativeness") |>
  mutate(VP    = factor(VP),
         model = factor(model, levels = c("PECAN", "nar_ex")))

max_score_red_pilot <- length(nw02_items_red_pilot) * 5

rep_desc_red_pilot <- rep_long_red_pilot |>
  group_by(model) |>
  summarise(
    M    = round(mean(representativeness), 2),
    SD   = round(sd(representativeness),   2),
    Min  = min(representativeness),
    Max  = max(representativeness),
    .groups = "drop"
  ) |>
  mutate(
    M_pct   = round(M   / max_score_red_pilot * 100, 1),
    SD_pct  = round(SD  / max_score_red_pilot * 100, 1),
    Min_pct = round(Min / max_score_red_pilot * 100, 1),
    Max_pct = round(Max / max_score_red_pilot * 100, 1)
  )

cat("\n[5b] Descriptives (reduced sum scores, max =", max_score_red_pilot, "| percentages of max):\n")
print(rep_desc_red_pilot)

# --- 5b.5 Assumption checks ---

cat("\n[5b] Shapiro-Wilk normality checks (reduced scores):\n")
for (cond in c("PECAN", "nar_ex")) {
  x  <- rep_complete_red_pilot[[cond]]
  sw <- shapiro.test(x)
  cat(sprintf("  %-8s W = %.3f, p = %.4f\n", cond, sw$statistic, sw$p.value))
}

# --- 5b.6 RM-ANOVA (reduced scores) ---
# With 2 conditions, sphericity is trivially met; RM-ANOVA = paired t-test.

rep_anova_red_pilot <- aov_ez(
  id          = "VP",
  dv          = "representativeness",
  data        = rep_long_red_pilot,
  within      = "model",
  anova_table = list(es = "pes")
)

cat("\n[5b] RM-ANOVA — representativeness, network-specific items only:\n")
print(rep_anova_red_pilot)

# --- 5b.7 Paired t-test (equivalent with 2 conditions; easier to report) ---

cat("\n[5b] Paired t-test (PECAN vs. nar_ex, reduced items):\n")
print(t.test(rep_complete_red_pilot$PECAN, rep_complete_red_pilot$nar_ex, paired = TRUE))

# --- 5b.7b Wilcoxon signed-rank (non-parametric sensitivity) ---

rep_wilcox_red_pilot <- wilcox.test(rep_complete_red_pilot$PECAN, rep_complete_red_pilot$nar_ex,
                               paired = TRUE, exact = FALSE)
cat("\n[5b] Wilcoxon signed-rank (sensitivity check):\n")
print(rep_wilcox_red_pilot)

wilcox_tbl_red_pilot <- tibble(
  pair          = "PECAN vs. nar_ex",
  V             = rep_wilcox_red_pilot$statistic,
  p_uncorrected = round(rep_wilcox_red_pilot$p.value, 4)
)
cat("\n"); print(wilcox_tbl_red_pilot)

# --- 5b.7c Profile plot (reduced scores) ---

ggplot(rep_long_red_pilot, aes(x = model, y = representativeness, group = VP)) +
  geom_line(colour = "grey70", linewidth = 0.4) +
  geom_point(colour = "grey50", size = 2) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               colour = "#2E86AB", linewidth = 1.3, linetype = "dashed") +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               shape = 18, size = 5, colour = "#2E86AB") +
  scale_x_discrete(labels = c("PECAN", "Nar-Ex")) +
  scale_y_continuous(limits = c(0, max_score_red_pilot), breaks = seq(0, max_score_red_pilot, 5)) +
  labs(
    x = "Network model",
    y = "Representativeness and feasibility (sum score)"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 13),
        axis.text    = element_text(size = 14))

ggsave("output/pilotstudy/plots/representativeness_profile_plot_reduced.png",
       width = 7, height = 5, dpi = 150)

# --- 5b.8 Save outputs ---

rep_posthoc_red_tbl_pilot <- as.data.frame(
  summary(pairs(emmeans::emmeans(rep_anova_red_pilot, ~ model), adjust = "bonferroni"))
) |> rename(pair = contrast, t_ratio = t.ratio, p_bonf = p.value)

write.csv(rep_data_red_pilot,
          "output/pilotstudy/tables/5b_representativeness_sumscores_reduced.csv",    row.names = FALSE)
write.csv(rep_desc_red_pilot,
          "output/pilotstudy/tables/5b_representativeness_descriptives_reduced.csv", row.names = FALSE)
write.csv(as.data.frame(rep_anova_red_pilot$anova_table),
          "output/pilotstudy/tables/5b_representativeness_anova_reduced.csv",        row.names = TRUE)
write.csv(rep_posthoc_red_tbl_pilot,
          "output/pilotstudy/tables/5b_representativeness_posthoc_reduced.csv",      row.names = FALSE)
write.csv(wilcox_tbl_red_pilot,
          "output/pilotstudy/tables/5b_representativeness_wilcox_pairs_reduced.csv", row.names = FALSE)

cat("Saved reduced-item representativeness tables (Section 5b).\n")

# --- Save all key results to a single RData file ---

save(
  t1_pilot, dass_data_pilot, alpha_result_pilot, alpha_result_t2_pilot, dass_t2_items_data_pilot,
  age_stats_pilot, gender_table_pilot, education_table_pilot,
  irr_results_pilot,
  jaccard_results_pilot, jaccard_summary_pilot,
  # Unregulated PECAN results
  edge_cor_results_pilot, centrality_cor_results_pilot,
  eif_pecan_pilot, edge_comparison_pilot,
  # Regulated PECAN results
  edge_cor_results_reg_pilot, edge_cor_combined_pilot,
  centrality_cor_results_reg_pilot, centrality_cor_combined_pilot,
  eif_pecan_reg_pilot, edge_comparison_reg_pilot,
  # Shared across both versions
  sif_pecan_pilot, sif_nar_pilot, eif_nar_pilot,
  loop_results_pilot, symptom_table_pilot,
  # Section 3c network metrics RM-ANOVA
  metric_data_long_pilot, anova_tbl_pilot, posthoc_tbl_pilot,
  # Section 5: representativeness (full item set)
  t2_pilot, dass_t2_data_pilot, rep_data_pilot, rep_desc_pilot, rep_anova_pilot,
  # Section 5b: representativeness (network-specific items only)
  rep_data_red_pilot, rep_desc_red_pilot, rep_anova_red_pilot, wilcox_tbl_red_pilot,
  file = "output/pilotstudy/N2N_results.RData"
)

# ===========================================================================
# SECTION 6: REPORTING OVERVIEW TABLES
# ===========================================================================

cat("\n=== SECTION 7: REPORTING OVERVIEW TABLES ===\n")

sumstat <- function(x, label, digits = 2) {
  x <- x[!is.na(x)]
  tibble(label = label, n = length(x),
         M   = round(mean(x), digits), SD  = round(sd(x),  digits),
         Min = round(min(x),  digits), Max = round(max(x), digits))
}

# ── Table 1: Sample descriptives (DASS-21) ───────────────────────────────────
tbl_sample_pilot <- bind_rows(
  sumstat(dass_data_pilot$DASS_total,       "DASS-21 total score — T1 (0–63)"),
  sumstat(dass_t2_data_pilot$DASS_T2_total, "DASS-21 total score — T2 (0–63)")
)
cat("\n[Table 1] Sample descriptives:\n"); print(tbl_sample_pilot)
write.csv(tbl_sample_pilot, "output/pilotstudy/tables/6_verview_01_sample_descriptives.csv",
          row.names = FALSE)

# ── Table 2: Network descriptives ────────────────────────────────────────────
# Node/edge counts, density, mean node severity, and mean edge width per
# network version. nar has no edge weights; strength = edge count for nar.
# Mean edge width (1–5) is PECAN-specific.

pilot_pecan_n_nodes        <- map_int(pecan_nodes_list_pilot,     nrow)
pilot_pecan_n_nodes_reg    <- map_int(pecan_nodes_list_pilot,     nrow)  # same: reg keeps all nodes
pilot_nar_n_nodes          <- map_int(nar_nodes_list_pilot,        nrow)
pilot_pecan_sev_per_vp     <- map_dbl(pecan_nodes_list_pilot,     ~ mean(.x$value, na.rm = TRUE))
pilot_pecan_mean_width     <- map_dbl(pecan_edges_list_pilot,     ~ mean(.x$width, na.rm = TRUE))
pilot_pecan_mean_width_reg <- map_dbl(pecan_edges_reg_list_pilot, ~ mean(.x$width, na.rm = TRUE))

# Compute edge counts, strength, and density inline (Section 4 removed)
.pilot_dens <- function(n_e, n_n) ifelse(n_n > 1, n_e / (n_n * (n_n - 1)), NA_real_)

pilot_pecan_n_edges_v   <- map_int(vps_pilot, ~ nrow(pecan_edges_list_pilot[[.x]]))
pilot_pecan_strength_v  <- map_dbl(vps_pilot, ~ sum(pecan_edges_list_pilot[[.x]]$width, na.rm = TRUE))
pilot_pecan_density_v   <- .pilot_dens(pilot_pecan_n_edges_v, pilot_pecan_n_nodes)
pilot_pecan_n_edges_reg_v  <- map_int(vps_pilot, ~ nrow(pecan_edges_reg_list_pilot[[.x]]))
pilot_pecan_strength_reg_v <- map_dbl(vps_pilot, ~ sum(pecan_edges_reg_list_pilot[[.x]]$width, na.rm = TRUE))
pilot_pecan_density_reg_v  <- .pilot_dens(pilot_pecan_n_edges_reg_v, pilot_pecan_n_nodes_reg)
pilot_nar_n_edges_v     <- map_int(vps_pilot, ~ nrow(distinct(nar_edges_list_pilot[[.x]], from, to)))
pilot_nar_density_v     <- .pilot_dens(pilot_nar_n_edges_v, pilot_nar_n_nodes)

tbl_networks_pilot <- bind_rows(
  # PECAN unregulated
  sumstat(pilot_pecan_n_nodes,           "PECAN (unreg) — nodes per participant"),
  sumstat(pilot_pecan_n_edges_v,         "PECAN (unreg) — edges per participant"),
  sumstat(pilot_pecan_strength_v,        "PECAN (unreg) — total edge weight (strength)"),
  sumstat(pilot_pecan_mean_width,        "PECAN (unreg) — mean edge width (1–5)"),
  sumstat(pilot_pecan_density_v,         "PECAN (unreg) — density"),
  sumstat(pilot_pecan_sev_per_vp,        "PECAN (unreg) — mean node severity (1–10)"),
  # PECAN regulated
  sumstat(pilot_pecan_n_nodes_reg,       "PECAN (reg)   — nodes per participant"),
  sumstat(pilot_pecan_n_edges_reg_v,     "PECAN (reg)   — edges per participant"),
  sumstat(pilot_pecan_strength_reg_v,    "PECAN (reg)   — total edge weight (strength)"),
  sumstat(pilot_pecan_mean_width_reg,    "PECAN (reg)   — mean edge width (1–5)"),
  sumstat(pilot_pecan_density_reg_v,     "PECAN (reg)   — density"),
  sumstat(pilot_pecan_sev_per_vp,        "PECAN (reg)   — mean node severity (1–10)"),
  # nar (full narrative network)
  sumstat(pilot_nar_n_nodes,             "nar — nodes per participant"),
  sumstat(pilot_nar_n_edges_v,           "nar — edges per participant"),
  sumstat(pilot_nar_density_v,           "nar — density"),
  # Feedback loops (cycles ≤ 4)
  sumstat(loop_results_pilot$pecan_le4,              "PECAN (unreg) — feedback loops (≤4) per participant"),
  sumstat(loop_results_pilot$pecan_reg_le4,          "PECAN (reg)   — feedback loops (≤4) per participant"),
  sumstat(loop_results_pilot$nar_le4,                "nar           — feedback loops (≤4) per participant")
)
cat("\n[Table 2] Network descriptives:\n"); print(tbl_networks_pilot, n = Inf)
write.csv(tbl_networks_pilot, "output/pilotstudy/tables/6_overview_02_network_descriptives.csv",
          row.names = FALSE)

cat("\n=== Analysis complete. All results saved to output/ ===\n")
cat("Tables: output/pilotstudy/tables/\n")
cat("Plots:  output/pilotstudy/plots/\n")
cat("RData:  output/pilotstudy/N2N_results.RData\n")