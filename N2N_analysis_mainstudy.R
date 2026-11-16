# =============================================================================
# N2N Main Study — Master's Thesis Data Analysis
# Comparing Narrative-to-Network (nar) with PECAN (gold standard)
# N = 50 participants (VP23–VP72)
# =============================================================================

# ===========================================================================
# SECTION 0: SETUP & DATA LOADING
# ===========================================================================

rm(list = ls())   # clear environment so all objects are rebuilt from source files

library(tidyverse)
library(psych)
library(qgraph)
library(PECAN2)
library(visNetwork)
library(irr)
library(emmeans)   # post-hoc contrasts for Section 5 RM-ANOVA
library(ppcor)     # partial Spearman correlations (combined script)
suppressWarnings(library(afex))

# ppcor loads MASS, which masks dplyr::select and dplyr::filter.
# Re-assert the dplyr versions explicitly.
select <- dplyr::select
filter <- dplyr::filter

# Restore default warning behaviour (warn = 1: print warnings, do not error).
# This prevents warn = 2 sessions from converting warnings to hard errors.
options(warn = 1)

# --- 0.1 Create working directories and copy raw data to data/bearbeitet/ ---
# IMPORTANT: All analysis reads from data/bearbeitet/. Raw files are never modified.

dir.create("data/bearbeitet/mainstudy/nar",        recursive = TRUE, showWarnings = FALSE)
dir.create("data/bearbeitet/mainstudy/nar_r2",     recursive = TRUE, showWarnings = FALSE)
dir.create("data/bearbeitet/mainstudy/PECAN",      recursive = TRUE, showWarnings = FALSE)
dir.create("output/mainstudy/tables",              recursive = TRUE, showWarnings = FALSE)
dir.create("output/mainstudy/plots",               recursive = TRUE, showWarnings = FALSE)
dir.create("output/mainstudy/networks",            recursive = TRUE, showWarnings = FALSE)

file.copy(
  "data/rawdata/mainstudy/N2N T1 Data mainstudy.csv",
  "data/bearbeitet/mainstudy/N2N T1 Data mainstudy.csv",
  overwrite = TRUE
)

nar_raw    <- list.files("data/rawdata/mainstudy/nar",    pattern = "\\.csv$", full.names = TRUE)
pecan_raw  <- list.files("data/rawdata/mainstudy/PECAN",  pattern = "\\.csv$", full.names = TRUE)
nar_r2_raw <- list.files("data/rawdata/mainstudy/nar_r2", pattern = "\\.csv$", full.names = TRUE)
file.copy(nar_raw,    "data/bearbeitet/mainstudy/nar/",    overwrite = TRUE)
file.copy(pecan_raw,  "data/bearbeitet/mainstudy/PECAN/",  overwrite = TRUE)
file.copy(nar_r2_raw, "data/bearbeitet/mainstudy/nar_r2/", overwrite = TRUE)
file.copy(
  "data/rawdata/mainstudy/N2N T2 Data mainstudy.csv",
  "data/bearbeitet/mainstudy/N2N T2 Data mainstudy.csv",
  overwrite = TRUE
)

# --- 0.2 Load T1 survey data ---
# File is UTF-16LE encoded and tab-separated (exported from SoSciSurvey)

t1 <- read.delim(
  "data/bearbeitet/mainstudy/N2N T1 Data mainstudy.csv",
  fileEncoding = "UTF-16LE",
  sep          = "\t"
)

# --- 0.3 Load per-participant network data into named lists ---

vps <- sprintf("VP%02d", 23:72)

nar_nodes_list   <- setNames(
  lapply(vps, \(vp) read.csv(paste0("data/bearbeitet/mainstudy/nar/",   vp, "_nar_nodes.csv"))),
  vps
)
nar_edges_list   <- setNames(
  lapply(vps, \(vp) read.csv(paste0("data/bearbeitet/mainstudy/nar/",   vp, "_nar_edges.csv"))),
  vps
)
pecan_nodes_list <- setNames(
  lapply(vps, \(vp) {
    df <- read.csv(paste0("data/bearbeitet/mainstudy/PECAN/", vp, "_PECAN_nodes.csv"))
    df$label <- trimws(df$label)   # remove leading/trailing whitespace (e.g. VP57 node 19)
    df
  }),
  vps
)
pecan_edges_list <- setNames(
  lapply(vps, \(vp) read.csv(paste0("data/bearbeitet/mainstudy/PECAN/", vp, "_PECAN_edges.csv"))),
  vps
)

# r2 data: one folder per rater, auto-detected from rawdata/mainstudy/nar_r2.*/
# Each folder is named e.g. nar_r2.1_G and contains *_nar_nodes_r2.csv / *_nar_edges_r2.csv.
# r2_raters: named list (folder name → list(vps, nodes, edges)); empty rater folders are skipped.
r2_rater_dirs <- list.dirs("data/rawdata/mainstudy", recursive = FALSE, full.names = FALSE)
r2_rater_dirs <- sort(r2_rater_dirs[grepl("^nar_r2\\.", r2_rater_dirs)])

r2_raters <- setNames(lapply(r2_rater_dirs, function(rd) {
  folder <- file.path("data/rawdata/mainstudy", rd)
  vps    <- list.files(folder, pattern = "_nar_nodes_r2\\.csv$") |>
    sub("_nar_nodes_r2\\.csv$", "", x = _)
  if (length(vps) == 0) return(NULL)
  list(
    vps   = vps,
    nodes = setNames(lapply(vps, \(vp) read.csv(file.path(folder, paste0(vp, "_nar_nodes_r2.csv")))), vps),
    edges = setNames(lapply(vps, \(vp) read.csv(file.path(folder, paste0(vp, "_nar_edges_r2.csv")))), vps)
  )
}), r2_rater_dirs)
r2_raters <- Filter(Negate(is.null), r2_raters)

cat("r2 rater folders found (n =", length(r2_raters), "):\n")
for (rd in names(r2_raters)) {
  cat("  ", rd, ": n =", length(r2_raters[[rd]]$vps), "VPs —",
      paste(r2_raters[[rd]]$vps, collapse = ", "), "\n")
}

n_vps <- length(vps)

# --- 0.3b Explicit subset ---
# Main study nar files contain an ex.im column ("e" = explicit, "i" = implicit).
# Nodes and edges tagged "e" form the explicit network (nar_ex).
# IMPORTANT: these lists must be built before the pooled objects below.

nar_ex_nodes_list <- lapply(nar_nodes_list, \(df) df |> filter(ex.im == "e"))
nar_ex_edges_list <- lapply(nar_edges_list, \(df) df |> filter(ex.im == "e"))

# Pool all nodes and edges across participants (used in Section 5)
all_nar_nodes   <- map_dfr(vps, ~ nar_nodes_list[[.x]]   |> mutate(vp = .x))
all_nar_edges   <- map_dfr(vps, ~ nar_edges_list[[.x]]   |> mutate(vp = .x))
all_pecan_nodes <- map_dfr(vps, ~ pecan_nodes_list[[.x]] |> mutate(vp = .x))
all_pecan_edges <- map_dfr(vps, ~ pecan_edges_list[[.x]] |> mutate(vp = .x))

# Pooled explicit subset (used in Section 4)
all_nar_ex_nodes <- map_dfr(vps, ~ nar_ex_nodes_list[[.x]] |> mutate(vp = .x))
all_nar_ex_edges <- map_dfr(vps, ~ nar_ex_edges_list[[.x]] |> mutate(vp = .x))

cat("Explicit node count (across all VPs):",
    nrow(map_dfr(nar_ex_nodes_list, \(x) x)), "\n")

# --- 0.4 PECAN edge regulation ---
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

pecan_edges_reg_list <- setNames(
  lapply(vps, \(vp) regulate_pecan_edges(pecan_edges_list[[vp]], pecan_nodes_list[[vp]])),
  vps
)
all_pecan_edges_reg <- map_dfr(vps, ~ pecan_edges_reg_list[[.x]] |> mutate(vp = .x))

# --- 0.5 Helper functions ---

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
#output from make_adj ist adj_nar_list and adj_pecan_list 
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


# Compare two per-participant Spearman rho vectors (paired test).
# Always uses Wilcoxon signed-rank (n is too small for Shapiro-Wilk to be informative).
# Returns a one-row tibble.
compare_rho_fn <- function(rho_a, rho_b, label) {
  valid <- !is.na(rho_a) & !is.na(rho_b)
  a_v   <- rho_a[valid]; b_v <- rho_b[valid]; d <- a_v - b_v
  n     <- length(d)
  if (n < 3)
    return(tibble(comparison = label, test = NA_character_,
                  statistic = NA_real_, df = NA_real_,
                  p_value = NA_real_, n = n,
                  mean_diff = NA_real_))
  wt <- wilcox.test(a_v, b_v, paired = TRUE, exact = FALSE)
  tibble(comparison = label, test = "Wilcoxon signed-rank",
         statistic  = round(wt$statistic, 3), df = NA_real_,
         p_value    = round(wt$p.value,   4), n = n,
         mean_diff  = round(mean(d),      3))
}


# Pre-compute co-occurrence: for each VP, which directed node pairs co-occur?
#all possible edges for each VP are listed --> used to compare the actual edges later on
build_cooccurrence <- function(nodes_list, vps) {
  map_dfr(vps, function(vp) {
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
age_stats <- t1 |>
  summarise(
    mean_age = round(mean(SD02_01, na.rm = TRUE), 2),
    sd_age   = round(sd(SD02_01,   na.rm = TRUE), 2),
    min_age  = min(SD02_01, na.rm = TRUE),
    max_age  = max(SD02_01, na.rm = TRUE)
  )
cat("Age — Mean:", age_stats$mean_age, "| SD:", age_stats$sd_age,
    "| Range:", age_stats$min_age, "–", age_stats$max_age, "\n")

# Gender (SD05): 1 = female, 2 = male
gender_table <- t1 |>
  mutate(gender = case_when(
    SD05 == 1 ~ "female",
    SD05 == 2 ~ "male",
    TRUE      ~ "other/missing"
  )) |>
  count(gender) |>
  mutate(pct = round(n / sum(n) * 100, 1))
cat("\nGender distribution:\n")
print(gender_table)

# Education (SD11)
education_table <- t1 |>
  count(SD11) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  arrange(SD11)
cat("\nEducation distribution (SD11):\n")
print(education_table)

write.csv(as.data.frame(age_stats), "output/mainstudy/tables/1_descriptives_age.csv",       row.names = FALSE)
write.csv(gender_table,             "output/mainstudy/tables/1_descriptives_gender.csv",    row.names = FALSE)
write.csv(education_table,          "output/mainstudy/tables/1_descriptives_education.csv", row.names = FALSE)


# ===========================================================================
# SECTION 2: RELIABILITY & CONCORDANCE
# ===========================================================================
#relevant for the master thesis were the following output tables:
#2b_cronbach_alpha_DASS21 #internal consistency of the DASS-21
#2c_irr_global_kappa_nar_r2.1_G #interrater-reliability measured between main rater and second rater 1 with a reference set pooled across all participants
#2c_irr_global_kappa_nar_r2.2_A #interrater-reliability measured between main rater and second rater 2 with a reference set pooled across all participants
#2c_irr_global_kappa_nar_r2.3_F #interrater-reliability measured between main rater and second rater 3 with a reference set pooled across all participants
#2c_irr_global_kappa_nar_r2.4_J #interrater-reliability measured between main rater and second rater 4 with a reference set pooled across all participants
#2d_cronbach_alpha_representativeness #internal consistency of the representativeness and feasibility scale

cat("\n=== SECTION 2: RELIABILITY & CONCORDANCE ===\n")

# ---------------------------------------------------------------------------
# --- 2a. DASS-21 sum score (recode 1–4 → 0–3 per standard scoring) ---
# ---------------------------------------------------------------------------

dass_items <- paste0("DA01_", sprintf("%02d", 1:21))

dass_data <- t1 |>
  select(VP01_01, all_of(dass_items)) |>
  # Recode: original 1–4 → standard 0–3
  mutate(across(all_of(dass_items), ~ .x - 1)) |>
  mutate(DASS_total = rowSums(across(all_of(dass_items)), na.rm = FALSE))

cat("DASS-21 sum scores per participant:\n")
print(dass_data |> select(VP01_01, DASS_total))
cat("DASS-21 total — Mean:", round(mean(dass_data$DASS_total, na.rm = TRUE), 2),
    "| SD:", round(sd(dass_data$DASS_total, na.rm = TRUE), 2),
    "| Range:", min(dass_data$DASS_total, na.rm = TRUE), "–",
    max(dass_data$DASS_total, na.rm = TRUE), "\n")

write.csv(dass_data |> select(VP01_01, DASS_total),
          "output/mainstudy/tables/2a_DASS21_scores.csv", row.names = FALSE)

# --- 2a-T2. DASS-21 T2 sum scores ---
# T2 file was copied to bearbeitet/ in Section 0.1
# Main study T2 is semicolon-separated (SoSciSurvey export format)

t2 <- read.delim(
  "data/bearbeitet/mainstudy/N2N T2 Data mainstudy.csv",
  fileEncoding = "UTF-16LE",
  sep = ";"
)

dass_t2_items <- paste0("DA01_", sprintf("%02d", 1:21))
dass_t2_data  <- t2 |>
  select(vp = VP, all_of(dass_t2_items)) |>
  mutate(across(all_of(dass_t2_items), ~ .x - 1)) |>
  mutate(DASS_T2_total = rowSums(across(all_of(dass_t2_items)), na.rm = FALSE)) |>
  select(vp, DASS_T2_total)

cat("\nDASS-21 T2 sum scores:\n")
print(dass_t2_data)
cat("T2 DASS — Mean:", round(mean(dass_t2_data$DASS_T2_total, na.rm = TRUE), 2),
    "| SD:", round(sd(dass_t2_data$DASS_T2_total, na.rm = TRUE), 2),
    "| Missing:", sum(is.na(dass_t2_data$DASS_T2_total)), "\n")

write.csv(dass_t2_data, "output/mainstudy/tables/2a_DASS21_T2_scores.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 2b. Cronbach's alpha for DASS-21 (T1 and T2) ---
# ---------------------------------------------------------------------------

# T1: all 50 participants
alpha_result <- withCallingHandlers(
  psych::alpha(dass_data |> select(all_of(dass_items)), check.keys = TRUE),
  warning = function(w) invokeRestart("muffleWarning")
)
cat("\nCronbach's Alpha for DASS-21 — T1 (n = 50):\n")
print(alpha_result$total)

# T2: only participants with complete DASS item responses
dass_t2_items_data <- t2 |>
  select(all_of(dass_t2_items)) |>
  mutate(across(everything(), ~ .x - 1)) |>
  filter(complete.cases(pick(everything())))

alpha_result_t2 <- withCallingHandlers(
  psych::alpha(dass_t2_items_data, check.keys = TRUE),
  warning = function(w) invokeRestart("muffleWarning")
)
cat(sprintf("\nCronbach's Alpha for DASS-21 — T2 (n = %d):\n", nrow(dass_t2_items_data)))
print(alpha_result_t2$total)

# Save both
alpha_dass_combined <- bind_rows(
  as.data.frame(alpha_result$total)    |> mutate(timepoint = "T1", n = nrow(dass_data),         .before = 1),
  as.data.frame(alpha_result_t2$total) |> mutate(timepoint = "T2", n = nrow(dass_t2_items_data), .before = 1)
)
write.csv(alpha_dass_combined, "output/mainstudy/tables/2b_cronbach_alpha_DASS21.csv",
          row.names = FALSE)

# --- 2c. Inter-rater Reliability (nar r1 vs. nar r2) — looped per rater ---
#
# r2_raters is a named list built in Section 0 from rawdata/mainstudy/nar_r2.*/ folders.
# Each rater is processed independently so that IRR estimates are not conflated across raters.
#
# Per rater, two sets of metrics are computed:
#   (i)  Nodes & edges (binary presence/absence) → Cohen's Kappa (irr::kappa2)
#         + Jaccard index (node sets)
#         + ICC two-way consistency & Cronbach's alpha for in-/out-degree
#   (ii) Kappa with a within-rater global reference set
#         (node/edge vocab = union across all VPs of that rater)
#
# Per-rater CSVs: irr_<rater_id>.csv, irr_global_kappa_<rater_id>.csv
# Combined CSVs:  irr_all_raters.csv, irr_global_kappa_all_raters.csv

if (length(r2_raters) == 0) {
  cat("\nSection 2c skipped: no r2 rater data available for main study.\n")
} else {

all_irr          <- list()
all_global       <- list()
deg_fn_list      <- list()
build_adj_list   <- list()
global_ref_list  <- list()

for (rater_id in names(r2_raters)) {

  vps_r2            <- r2_raters[[rater_id]]$vps
  nar_nodes_r2_list <- r2_raters[[rater_id]]$nodes
  nar_edges_r2_list <- r2_raters[[rater_id]]$edges

  cat(sprintf("\n=== IRR: Rater %s (n = %d VPs: %s) ===\n",
              rater_id, length(vps_r2), paste(vps_r2, collapse = ", ")))

  # --- (i) Main: kappa + Jaccard + ICC/alpha ---
  irr_results <- map_dfr(vps_r2, function(vp) {
    n1 <- nar_nodes_list[[vp]]
    e1 <- distinct(nar_edges_list[[vp]], from, to)
    n2 <- nar_nodes_r2_list[[vp]]
    e2 <- distinct(nar_edges_r2_list[[vp]], from, to)

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
    if (length(shared_ids) >= 3) {
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

  cat("\nIRR (kappa / ICC / Jaccard) —", rater_id, ":\n")
  print(irr_results)
  cat("\n--- Summary (mean ± SD | Range) ---\n")
  cat("Nodes  — Kappa:       ", fmt(irr_results$kappa_nodes),   "| Range:", rng(irr_results$kappa_nodes),   "\n")
  cat("Edges  — Kappa:       ", fmt(irr_results$kappa_edges),   "| Range:", rng(irr_results$kappa_edges),   "\n")
  cat("Nodes  — Jaccard:     ", fmt(irr_results$jaccard_nodes), "| Range:", rng(irr_results$jaccard_nodes), "\n")
  cat("In-deg — ICC:         ", fmt(irr_results$icc_indeg),     "| Range:", rng(irr_results$icc_indeg),     "\n")
  cat("In-deg — Cronbach α:  ", fmt(irr_results$alpha_indeg),   "| Range:", rng(irr_results$alpha_indeg),   "\n")
  cat("Out-deg — ICC:        ", fmt(irr_results$icc_outdeg),    "| Range:", rng(irr_results$icc_outdeg),    "\n")
  cat("Out-deg — Cronbach α: ", fmt(irr_results$alpha_outdeg),  "| Range:", rng(irr_results$alpha_outdeg),  "\n")
  write.csv(
    append_summary(irr_results,
      c("kappa_nodes", "kappa_edges", "jaccard_nodes",
        "icc_indeg", "alpha_indeg", "icc_outdeg", "alpha_outdeg")),
    sprintf("output/mainstudy/tables/2c_irr_%s.csv", rater_id), row.names = FALSE
  )

  # Save build_adj output per VP for inspection
  build_adj_list[[rater_id]] <- setNames(
    lapply(vps_r2, function(vp) {
      e1      <- distinct(nar_edges_list[[vp]], from, to)
      e2      <- distinct(nar_edges_r2_list[[vp]], from, to)
      all_ids <- sort(union(
        as.character(nar_nodes_list[[vp]]$id),
        as.character(nar_nodes_r2_list[[vp]]$id)
      ))
      list(
        rater1 = build_adj(e1, all_ids),
        rater2 = build_adj(e2, all_ids)
      )
    }),
    vps_r2
  )

  # Save deg_fn output per VP for inspection
  deg_fn_list[[rater_id]] <- setNames(
    lapply(vps_r2, function(vp) {
      e1         <- distinct(nar_edges_list[[vp]], from, to)
      e2         <- distinct(nar_edges_r2_list[[vp]], from, to)
      shared_ids <- intersect(
        as.character(nar_nodes_list[[vp]]$id),
        as.character(nar_nodes_r2_list[[vp]]$id)
      )
      if (length(shared_ids) < 3) return(NULL)
      list(
        rater1 = deg_fn(e1, shared_ids),
        rater2 = deg_fn(e2, shared_ids)
      )
    }),
    vps_r2
  )

  # --- (ii) Kappa with within-rater global reference set ---
  # Global vocab = union across all VPs rated by this rater.
  # Makes kappa values structurally comparable within the rater's VP set.
  global_node_ids_r <- sort(unique(c(
    unlist(lapply(vps_r2, \(vp) as.character(nar_nodes_list[[vp]]$id))),
    unlist(lapply(vps_r2, \(vp) as.character(nar_nodes_r2_list[[vp]]$id)))
  )))
  global_edge_ids_r <- sort(unique(c(
    unlist(lapply(vps_r2, \(vp) {
      e <- distinct(nar_edges_list[[vp]], from, to)
      paste(e$from, e$to, sep = "->")
    })),
    unlist(lapply(vps_r2, \(vp) {
      e <- distinct(nar_edges_r2_list[[vp]], from, to)
      paste(e$from, e$to, sep = "->")
    }))
  )))

  # Save reference set for inspection
  global_ref_list[[rater_id]] <- list(
    node_ids = global_node_ids_r,
    edge_ids = global_edge_ids_r,
    vectors  = setNames(
      lapply(vps_r2, function(vp) {
        ids1 <- as.character(nar_nodes_list[[vp]]$id)
        ids2 <- as.character(nar_nodes_r2_list[[vp]]$id)
        e1   <- distinct(nar_edges_list[[vp]], from, to)
        e2   <- distinct(nar_edges_r2_list[[vp]], from, to)
        ep1  <- paste(e1$from, e1$to, sep = "->")
        ep2  <- paste(e2$from, e2$to, sep = "->")
        list(
          gnode_r1 = as.integer(global_node_ids_r %in% ids1),
          gnode_r2 = as.integer(global_node_ids_r %in% ids2),
          gedge_r1 = as.integer(global_edge_ids_r %in% ep1),
          gedge_r2 = as.integer(global_edge_ids_r %in% ep2)
        )
      }),
      vps_r2
    )
  )

  irr_global_kappa <- map_dfr(vps_r2, function(vp) {
    ids1 <- as.character(nar_nodes_list[[vp]]$id)
    ids2 <- as.character(nar_nodes_r2_list[[vp]]$id)
    e1   <- distinct(nar_edges_list[[vp]], from, to)
    e2   <- distinct(nar_edges_r2_list[[vp]], from, to)
    ep1  <- paste(e1$from, e1$to, sep = "->")
    ep2  <- paste(e2$from, e2$to, sep = "->")

    gnode_r1 <- as.integer(global_node_ids_r %in% ids1)
    gnode_r2 <- as.integer(global_node_ids_r %in% ids2)
    kappa_nodes_global <- if (sd(gnode_r1) == 0 && sd(gnode_r2) == 0) NA_real_ else
      safe_irr(round(irr::kappa2(data.frame(r1 = gnode_r1, r2 = gnode_r2))$value, 3))

    gedge_r1 <- as.integer(global_edge_ids_r %in% ep1)
    gedge_r2 <- as.integer(global_edge_ids_r %in% ep2)
    kappa_edges_global <- if (sd(gedge_r1) == 0 && sd(gedge_r2) == 0) NA_real_ else
      safe_irr(round(irr::kappa2(data.frame(r1 = gedge_r1, r2 = gedge_r2))$value, 3))

    tibble(vp = vp, kappa_nodes_global = kappa_nodes_global, kappa_edges_global = kappa_edges_global)
  })

  cat("\nIRR (global ref kappa) —", rater_id, ":\n")
  print(irr_global_kappa)
  cat("Nodes — Kappa (global ref):", fmt(irr_global_kappa$kappa_nodes_global),
      "| Range:", rng(irr_global_kappa$kappa_nodes_global), "\n")
  cat("Edges — Kappa (global ref):", fmt(irr_global_kappa$kappa_edges_global),
      "| Range:", rng(irr_global_kappa$kappa_edges_global), "\n")
  cat("(Reference: n =", length(global_node_ids_r), "node IDs,",
      length(global_edge_ids_r), "directed edge pairs)\n")
  write.csv(
    append_summary(irr_global_kappa, c("kappa_nodes_global", "kappa_edges_global")),
    sprintf("output/mainstudy/tables/2c_irr_global_kappa_%s.csv", rater_id), row.names = FALSE
  )

  # Collect for cross-rater combined tables
  all_irr[[rater_id]]    <- irr_results      |> mutate(rater = rater_id, .before = 1)
  all_global[[rater_id]] <- irr_global_kappa |> mutate(rater = rater_id, .before = 1)

} # end for rater_id

# --- Cross-rater combined tables ---
irr_all_raters <- bind_rows(all_irr)
irr_all_global <- bind_rows(all_global)

cat("\n=== Combined IRR: mean per rater (kappa / ICC / Jaccard) ===\n")
irr_all_raters |>
  group_by(rater) |>
  summarise(
    n_vps       = n(),
    kappa_nodes = round(mean(kappa_nodes,   na.rm = TRUE), 3),
    kappa_edges = round(mean(kappa_edges,   na.rm = TRUE), 3),
    jaccard     = round(mean(jaccard_nodes, na.rm = TRUE), 3),
    icc_indeg   = round(mean(icc_indeg,     na.rm = TRUE), 3),
    icc_outdeg  = round(mean(icc_outdeg,    na.rm = TRUE), 3),
    .groups = "drop"
  ) |> print()

cat("\n=== Combined IRR: overall mean, SD, and range across all raters ===\n")
irr_summary_metrics <- c("kappa_nodes", "kappa_edges", "jaccard_nodes",
                          "icc_indeg",   "alpha_indeg",
                          "icc_outdeg",  "alpha_outdeg")
irr_overall_summary <- map_dfr(irr_summary_metrics, function(m) {
  x <- irr_all_raters[[m]][!is.na(irr_all_raters[[m]])]
  tibble(metric = m,
         n      = length(x),
         mean   = round(mean(x), 3),
         sd     = round(sd(x),   3),
         min    = round(min(x),  3),
         max    = round(max(x),  3))
})
print(irr_overall_summary)

write.csv(irr_overall_summary,
  "output/mainstudy/tables/2c_irr_overall_summary.csv", row.names = FALSE)
write.csv(irr_all_raters,
  "output/mainstudy/tables/2c_irr_all_raters.csv", row.names = FALSE)
write.csv(irr_all_global,
  "output/mainstudy/tables/2c_irr_global_kappa_all_raters.csv", row.names = FALSE)

} # end if (length(r2_raters) > 0)

# ---------------------------------------------------------------------------
# --- 2d. Internal consistency of representativeness scales (NW02, NW06, NW07) ---
# T2 contains three 14-item representativeness questionnaires:
#   NW02 = PECAN network,  NW06 = nar network (full: ex + im),  NW07 = nar_ex network
# Items NW02_14, NW02_15, NW06_14, NW06_15, NW07_14, NW07_15 are reverse-scored.
# Cronbach's alpha is computed for each scale over all T2 respondents with
# complete data (NAs excluded per scale).
# ---------------------------------------------------------------------------

nw_rev_items <- c("NW02_14", "NW02_15", "NW06_14", "NW06_15", "NW07_14", "NW07_15")

# Recode reverse items: scale is 1–5, so reversed = 6 - x
t2_scored <- t2 |>
  mutate(across(all_of(nw_rev_items), ~ 6L - .x))

alpha_nw <- lapply(
  list(PECAN = "NW02", nar = "NW06", nar_ex = "NW07"),
  function(prefix) {
    items <- names(t2_scored)[grepl(paste0("^", prefix, "_"), names(t2_scored))]
    dat   <- t2_scored |> select(all_of(items)) |> filter(complete.cases(pick(everything())))
    withCallingHandlers(
      psych::alpha(dat, check.keys = FALSE),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }
)

cat("\nCronbach's Alpha — representativeness scales:\n")
for (nm in names(alpha_nw)) {
  a <- round(alpha_nw[[nm]]$total$raw_alpha, 3)
  n <- length(alpha_nw[[nm]]$scores)
  cat(sprintf("  %-8s (NW%s): α = %.3f  (n = %d)\n",
    nm, switch(nm, PECAN = "02", nar = "06", nar_ex = "07"), a, n))
}

alpha_nw_tbl <- map_dfr(names(alpha_nw), \(nm)
  as.data.frame(alpha_nw[[nm]]$total) |> mutate(scale = nm, .before = 1)
)
write.csv(alpha_nw_tbl, "output/mainstudy/tables/2d_cronbach_alpha_representativeness.csv",
          row.names = FALSE)


# ===========================================================================
# SECTION 3: STRUCTURAL VALIDITY
# ===========================================================================
#relevant for the master thesis were the following output tables:
#3a_jaccard_summary_by_network #jaccard indices of node overlap between nar networks and PECAN network
#3a_jaccard_pairwise_tests #tests whether nar_full and nar_ex differ significantly in their overlap with PECAN network
#3b_edge_concordance_spearman_combined #summary of Spearman correlations between narrative networks and PECAN networks
#3b_edge_concordance_nar_vs_ex_comparison #tests whether nar_full and nar_ex differ significantly in their overlap with PECAN networks
#3c_centrality_spearman_combined.csv #summary of Spearman correlations between narrative networks and PECAN networks
#3c_centrality_nar_vs_ex_comparison.csv #tests whether nar_full and nar_ex differ significantly in their overlap with PECAN networks
#3c_top_node_agreement_overview #asks the question whether nar and PECAN networks identified the same nodes to have the highest in/outdegree
#3d_node_value_dich_concordance_per_participant #PECAN node values are dichotomized using a per-participant median split; #summary of Spearman correlations between narrative networks and PECAN networks
#3d_node_value_dich_nar_vs_ex_comparison #tests whether nar_full and nar_ex differ significantly in their overlap with PECAN network

cat("\n=== SECTION 3: STRUCTURAL VALIDITY ===\n")

# ---------------------------------------------------------------------------
# --- 3a. Jaccard Index (node concordance per participant) ---
# ---------------------------------------------------------------------------
# Computed for nar (full) and nar_ex (explicit only) vs. PECAN (unregulated
# and regulated). PECAN regulation only trims edges, never nodes, so the node
# sets — and therefore Jaccard values — are identical between regulated and
# unregulated PECAN. The primary contrast of interest is nar vs. nar_ex.
#
# Four conditions: nar × PECAN_unreg, nar × PECAN_reg,
#                  nar_ex × PECAN_unreg, nar_ex × PECAN_reg.
pecan_nodes_reg_list <- pecan_nodes_list   # regulation does not alter node sets

jaccard_fn <- function(nar_nodes_lst, nar_label, pecan_nodes_lst, pecan_label) {
  tibble(
    vp            = vps,
    nar_version   = nar_label,
    pecan_version = pecan_label,
    n_nar         = map_int(vps, ~ nrow(nar_nodes_lst[[.x]])),
    n_pecan       = map_int(vps, ~ nrow(pecan_nodes_lst[[.x]])),
    n_shared      = map_int(vps, ~ length(intersect(
      as.character(nar_nodes_lst[[.x]]$id),
      as.character(pecan_nodes_lst[[.x]]$id)
    ))),
    jaccard       = map_dbl(vps, ~ jaccard(
      as.character(nar_nodes_lst[[.x]]$id),
      as.character(pecan_nodes_lst[[.x]]$id)
    ))
  ) |>
    mutate(
      prop_shared_of_nar   = round(n_shared / n_nar,   3),
      prop_shared_of_pecan = round(n_shared / n_pecan, 3)
    )
}

jacc_nar_unreg <- jaccard_fn(nar_nodes_list,    "nar",    pecan_nodes_list,     "PECAN_unreg")
jacc_nar_reg   <- jaccard_fn(nar_nodes_list,    "nar",    pecan_nodes_reg_list, "PECAN_reg")
jacc_ex_unreg  <- jaccard_fn(nar_ex_nodes_list, "nar_ex", pecan_nodes_list,     "PECAN_unreg")
jacc_ex_reg    <- jaccard_fn(nar_ex_nodes_list, "nar_ex", pecan_nodes_reg_list, "PECAN_reg")

jaccard_results      <- bind_rows(jacc_nar_unreg, jacc_nar_reg, jacc_ex_unreg, jacc_ex_reg)
jaccard_reg_combined <- jaccard_results   # alias kept for downstream compatibility

jaccard_summary <- jaccard_results |>
  group_by(nar_version, pecan_version) |>
  summarise(
    mean_J               = round(mean(jaccard,               na.rm = TRUE), 3),
    sd_J                 = round(sd(jaccard,                 na.rm = TRUE), 3),
    min_J                = round(min(jaccard,                na.rm = TRUE), 3),
    max_J                = round(max(jaccard,                na.rm = TRUE), 3),
    mean_prop_of_nar     = round(mean(prop_shared_of_nar,    na.rm = TRUE), 3),
    mean_prop_of_pecan   = round(mean(prop_shared_of_pecan,  na.rm = TRUE), 3),
    .groups = "drop"
  )
jaccard_reg_summary <- jaccard_summary    # alias kept for downstream compatibility

cat("\nJaccard Index per participant (nar / nar_ex × PECAN unreg/reg):\n")
print(jaccard_results |> select(vp, nar_version, pecan_version, n_nar, n_pecan, n_shared, jaccard))
cat("\nJaccard summary by condition:\n")
print(jaccard_summary)

write.csv(jaccard_results, "output/mainstudy/tables/3a_jaccard_per_participant.csv",   row.names = FALSE)
write.csv(jaccard_summary, "output/mainstudy/tables/3a_jaccard_summary_by_network.csv", row.names = FALSE)

# --- 3a-CMP. Pairwise comparison of Jaccard conditions ---
# PECAN regulation does not alter node sets, so the only non-redundant
# comparison is nar vs. nar_ex. No Bonferroni correction needed (single test).
# Normality of paired differences checked with Shapiro-Wilk;
# paired t-test if p > .05, Wilcoxon signed-rank otherwise.

jacc_pairs <- list(
  list(a = jacc_nar_unreg$jaccard, b = jacc_ex_unreg$jaccard,
       label = "nar vs. nar_ex")
)

jacc_pair_tests <- map_dfr(jacc_pairs, function(pr) {
  valid <- !is.na(pr$a) & !is.na(pr$b)
  a_v   <- pr$a[valid]
  b_v   <- pr$b[valid]
  d     <- a_v - b_v
  n     <- length(d)
  if (n < 3 || sd(d) == 0) {
    return(tibble(
      comparison = pr$label, test = NA_character_,
      statistic = NA_real_, df = NA_real_,
      p_uncorr = NA_real_, p_bonf = NA_real_,
      n = n, mean_diff = round(mean(d), 3), sw_p = NA_real_
    ))
  }
  sw_p <- shapiro.test(d)$p.value
  if (sw_p > 0.05) {
    tt <- t.test(a_v, b_v, paired = TRUE)
    tibble(
      comparison = pr$label, test = "paired t-test",
      statistic  = round(tt$statistic,  3),
      df         = round(tt$parameter,  1),
      p_uncorr   = round(tt$p.value,    4),      p_bonf     = NA_real_,
      n          = n,
      mean_diff  = round(tt$estimate,   3),
      sw_p       = round(sw_p,          4)
    )
  } else {
    wt <- wilcox.test(a_v, b_v, paired = TRUE, exact = FALSE)
    tibble(
      comparison = pr$label, test = "Wilcoxon signed-rank",
      statistic  = round(wt$statistic,  3),
      df         = NA_real_,
      p_uncorr   = round(wt$p.value,    4),
      p_bonf     = NA_real_,
      n          = n,
      mean_diff  = round(mean(d),       3),
      sw_p       = round(sw_p,          4)
    )
  }
})

cat("\n[3a-CMP] nar vs. nar_ex Jaccard (single paired test):\n")
print(jacc_pair_tests)
write.csv(jacc_pair_tests, "output/mainstudy/tables/3a_jaccard_pairwise_tests.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Section 3 helper functions (used in 3b onward)
# ---------------------------------------------------------------------------

# Per-VP edge concordance: Spearman rho on binary adjacency matrices over
# shared nodes between a nar variant and a PECAN edge list.
edge_concord_fn <- function(nar_nodes_lst, nar_edges_lst, pecan_edges_lst,
                            list_name = NULL,
                            nar_mat_list_name = NULL, pecan_mat_list_name = NULL) {
  if (!is.null(list_name))           assign(list_name,           list(), envir = .GlobalEnv)
  if (!is.null(nar_mat_list_name))   assign(nar_mat_list_name,   list(), envir = .GlobalEnv)
  if (!is.null(pecan_mat_list_name)) assign(pecan_mat_list_name, list(), envir = .GlobalEnv)
  map_dfr(vps, function(vp) {
    shared_ids <- intersect(nar_nodes_lst[[vp]]$id, pecan_nodes_list[[vp]]$id)
    if (length(shared_ids) < 3)
      return(tibble(vp = vp, rho = NA_real_, p = NA_real_, n_shared_nodes = length(shared_ids)))
    adj_nar   <- make_adj(nar_edges_lst[[vp]], shared_ids)
    adj_pecan <- make_adj(pecan_edges_lst[[vp]], shared_ids)
    adj_pecan[adj_pecan > 0] <- 1
    if (!is.null(nar_mat_list_name)) {
      tmp <- get(nar_mat_list_name, envir = .GlobalEnv)
      tmp[[vp]] <- adj_nar
      assign(nar_mat_list_name, tmp, envir = .GlobalEnv)
    }
    if (!is.null(pecan_mat_list_name)) {
      tmp <- get(pecan_mat_list_name, envir = .GlobalEnv)
      tmp[[vp]] <- adj_pecan
      assign(pecan_mat_list_name, tmp, envir = .GlobalEnv)
    }
    idx       <- which(row(adj_nar) != col(adj_nar))
    nar_vec   <- adj_nar[idx]; pecan_vec <- adj_pecan[idx]
    if (!is.null(list_name)) {
      tmp <- get(list_name, envir = .GlobalEnv)
      tmp[[vp]] <- list(nar_vec = nar_vec, pecan_vec = pecan_vec)
      assign(list_name, tmp, envir = .GlobalEnv)
    }
    if (sd(nar_vec) == 0 || sd(pecan_vec) == 0)
      return(tibble(vp = vp, rho = NA_real_, p = NA_real_, n_shared_nodes = length(shared_ids)))
    ct <- cor.test(nar_vec, pecan_vec, method = "spearman", exact = FALSE)
    tibble(vp = vp, rho = round(ct$estimate, 3), p = round(ct$p.value, 4),
           n_shared_nodes = length(shared_ids))
  })
}

# Per-VP centrality Spearman: nar degree vs PECAN centrality over shared nodes.
centrality_fn <- function(nar_nodes_lst, nar_edges_lst, pecan_edges_lst, list_name = NULL) {
  if (!is.null(list_name)) assign(list_name, list(), envir = .GlobalEnv)
  map_dfr(vps, function(vp) {
    pn <- pecan_nodes_list[[vp]]; pe <- pecan_edges_lst[[vp]]
    ne <- nar_edges_lst[[vp]]
    na_row <- tibble(vp = vp,
                     rho_in_unwt=NA_real_, p_in_unwt=NA_real_, rho_in_wt=NA_real_,
                     rho_out_unwt=NA_real_, p_out_unwt=NA_real_, rho_out_wt=NA_real_,
                     rho_sevwod=NA_real_, n_shared=NA_integer_)
    pecan_cen <- tryCatch(
      pecanCen(pe, pn, centrality_by = "width", sev_weighted = "value"),
      error = function(e) NULL)
    if (is.null(pecan_cen)) return(na_row)
    pecan_cen <- pecan_cen |> mutate(id = as.character(id))
    nar_out    <- ne |> count(from, name = "out_deg") |> rename(id = from)
    nar_in     <- ne |> count(to,   name = "in_deg")  |> rename(id = to)
    nar_degree <- full_join(nar_out, nar_in, by = "id") |>
      mutate(across(c(out_deg, in_deg), ~ replace_na(.x, 0L)), id = as.character(id))
    if (!is.null(list_name)) {
      tmp <- get(list_name, envir = .GlobalEnv)
      tmp[[vp]] <- nar_degree
      assign(list_name, tmp, envir = .GlobalEnv)
    }
    shared <- inner_join(
      nar_degree |> select(id, in_deg, out_deg),
      pecan_cen  |> select(id, in_connect, out_connect, in_degree, out_degree, sevwod_combined),
      by = "id")
    if (nrow(shared) < 3) return(na_row |> mutate(n_shared = nrow(shared)))
    # Only use nodes that belong to nar_nodes_lst for this VP
    nar_ids <- as.character(nar_nodes_lst[[vp]]$id)
    shared  <- shared |> filter(id %in% nar_ids)
    if (nrow(shared) < 3) return(na_row |> mutate(n_shared = nrow(shared)))
    r_in_unwt  <- rho_fn(shared$in_deg,  shared$in_connect)
    r_in_wt    <- rho_fn(shared$in_deg,  shared$in_degree)
    r_out_unwt <- rho_fn(shared$out_deg, shared$out_connect)
    r_out_wt   <- rho_fn(shared$out_deg, shared$out_degree)
    r_sevwod   <- rho_fn(shared$out_deg, shared$sevwod_combined)
    tibble(vp=vp,
           rho_in_unwt=r_in_unwt$rho,  p_in_unwt=r_in_unwt$p,   rho_in_wt=r_in_wt$rho,
           rho_out_unwt=r_out_unwt$rho, p_out_unwt=r_out_unwt$p, rho_out_wt=r_out_wt$rho,
           rho_sevwod=r_sevwod$rho, n_shared=nrow(shared))
  })
}

# Print summary for an edge concordance result tibble
print_edge_concord <- function(res, label) {
  valid_rho <- res$rho[!is.na(res$rho)]
  cat(sprintf("\n%s — Edge concordance (Spearman rho, binary adj):\n", label))
  cat(sprintf("  n valid: %d / %d\n", length(valid_rho), nrow(res)))
  if (length(valid_rho) == 0) {
    cat("  One-sample Wilcoxon: not estimable (0 valid observations)\n")
    return(invisible(NULL))
  }
  cat("  Mean rho:", round(mean(valid_rho), 3),
      "| SD:", round(sd(valid_rho), 3),
      "| Range:", round(min(valid_rho), 3), "–",
      round(max(valid_rho), 3), "\n")
  cat("  One-sample Wilcoxon (H1: rho > 0): ")
  wt <- wilcox.test(valid_rho, mu = 0, alternative = "greater", exact = FALSE)
  cat("V =", wt$statistic, ", p =", round(wt$p.value, 4), "\n")
}

# Print summary for a centrality result tibble; returns a summary tibble invisibly.
print_centrality <- function(res, label) {
  rho_cen_cols_local <- c("rho_in_unwt","rho_in_wt","rho_out_unwt","rho_out_wt","rho_sevwod")
  cat(sprintf("\n%s — Centrality Spearman (nar degree vs. PECAN):\n", label))
  summary_rows <- map_dfr(rho_cen_cols_local, function(col) {
    valid_rho <- res[[col]][!is.na(res[[col]])]
    cat(sprintf("  %-40s Mean rho: %.3f | SD: %.3f | n = %d\n", col,
                mean(valid_rho), sd(valid_rho), length(valid_rho)))
    row <- tibble(
      label     = label,
      rho_col   = col,
      n_valid   = length(valid_rho),
      mean_rho  = round(mean(valid_rho), 3),
      sd_rho    = round(sd(valid_rho),   3),
      W         = NA_real_,
      p_value   = NA_real_
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

# ---------------------------------------------------------------------------
# --- 3b. Edge concordance (Spearman on binary adjacency matrices) ---
# For each VP: build binary adjacency matrices over concordant (shared) nodes,
# flatten to vectors, compute Spearman rho.
# Run for nar (full) and nar_ex vs. PECAN (unregulated and regulated).
# ---------------------------------------------------------------------------

edge_cor_results        <- edge_concord_fn(nar_nodes_list,    nar_edges_list,    pecan_edges_list,
                                           "edge_vec_list",        "adj_nar_list",       "adj_pecan_list")
edge_cor_results_ex     <- edge_concord_fn(nar_ex_nodes_list, nar_ex_edges_list, pecan_edges_list,
                                           "edge_vec_ex_list",     "adj_nar_ex_list",    "adj_pecan_ex_list")
edge_cor_results_reg    <- edge_concord_fn(nar_nodes_list,    nar_edges_list,    pecan_edges_reg_list,
                                           "edge_vec_reg_list",    "adj_nar_reg_list",   "adj_pecan_reg_list")
edge_cor_results_ex_reg <- edge_concord_fn(nar_ex_nodes_list, nar_ex_edges_list, pecan_edges_reg_list,
                                           "edge_vec_ex_reg_list", "adj_nar_ex_reg_list","adj_pecan_ex_reg_list")

print_edge_concord(edge_cor_results,        "nar (unregulated)")
print_edge_concord(edge_cor_results_ex,     "nar_ex (unregulated)")
print_edge_concord(edge_cor_results_reg,    "nar (regulated)")
print_edge_concord(edge_cor_results_ex_reg, "nar_ex (regulated)")

# Combined table: 2 nar variants × 2 PECAN versions
edge_cor_combined <- bind_rows(
  edge_cor_results        |> mutate(nar_version = "nar",    pecan_version = "unregulated"),
  edge_cor_results_ex     |> mutate(nar_version = "nar_ex", pecan_version = "unregulated"),
  edge_cor_results_reg    |> mutate(nar_version = "nar",    pecan_version = "regulated"),
  edge_cor_results_ex_reg |> mutate(nar_version = "nar_ex", pecan_version = "regulated")
) |> select(nar_version, pecan_version, everything())

write.csv(
  append_summary(edge_cor_combined, "rho", group_col = c("nar_version", "pecan_version")),
  "output/mainstudy/tables/3b_edge_concordance_spearman_combined.csv", row.names = FALSE
)
# Also save the original nar-only files for back-compatibility
write.csv(append_summary(edge_cor_results,     "rho"),
          "output/mainstudy/tables/3b_edge_concordance_spearman_unreg.csv",  row.names = FALSE)
write.csv(append_summary(edge_cor_results_reg, "rho"),
          "output/mainstudy/tables/3b_edge_concordance_spearman_reg.csv",    row.names = FALSE)

# --- 3b-CMP. nar vs. nar_ex edge concordance comparison ---
# Test whether nar (full, explicit + implicit) and nar_ex (explicit only) achieve
# significantly different per-VP Spearman rhos with PECAN.
# Two comparisons: unregulated and regulated PECAN.

edge_concord_cmp <- bind_rows(
  compare_rho_fn(edge_cor_results$rho,     edge_cor_results_ex$rho,
                 "nar vs. nar_ex — PECAN unregulated"),
  compare_rho_fn(edge_cor_results_reg$rho, edge_cor_results_ex_reg$rho,
                 "nar vs. nar_ex — PECAN regulated")
)

cat("\n--- 3a-CMP: nar vs. nar_ex edge concordance ---\n")
print(edge_concord_cmp)

write.csv(edge_concord_cmp,
          "output/mainstudy/tables/3b_edge_concordance_nar_vs_ex_comparison.csv",
          row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 3c. Node centrality comparison (Spearman) ---
# Per-participant Spearman rho between nar and PECAN node degree, computed
# separately for in-degree and out-degree (unweighted and weighted).
# Run for nar and nar_ex × unregulated/regulated PECAN.
#
# NOTE: Only participants with ≥ 3 shared nodes yield valid correlations;
# with n = 3–4 nodes per correlation, estimates are highly unstable.
# ---------------------------------------------------------------------------

rho_cen_cols <- c("rho_in_unwt","rho_in_wt","rho_out_unwt","rho_out_wt","rho_sevwod")

centrality_cor_results        <- centrality_fn(nar_nodes_list,    nar_edges_list,    pecan_edges_list,     "nar_degree_list")
centrality_cor_results_ex     <- centrality_fn(nar_ex_nodes_list, nar_ex_edges_list, pecan_edges_list,     "nar_degree_ex_list")
centrality_cor_results_reg    <- centrality_fn(nar_nodes_list,    nar_edges_list,    pecan_edges_reg_list, "nar_degree_reg_list")
centrality_cor_results_ex_reg <- centrality_fn(nar_ex_nodes_list, nar_ex_edges_list, pecan_edges_reg_list, "nar_degree_ex_reg_list")

centrality_wilcox_summary <- bind_rows(
  print_centrality(centrality_cor_results,        "nar (unregulated)"),
  print_centrality(centrality_cor_results_ex,     "nar_ex (unregulated)"),
  print_centrality(centrality_cor_results_reg,    "nar (regulated)"),
  print_centrality(centrality_cor_results_ex_reg, "nar_ex (regulated)")
)

# Combined table: 2 nar variants × 2 PECAN versions
centrality_cor_combined <- bind_rows(
  centrality_cor_results        |> mutate(nar_version="nar",    pecan_version="unregulated"),
  centrality_cor_results_ex     |> mutate(nar_version="nar_ex", pecan_version="unregulated"),
  centrality_cor_results_reg    |> mutate(nar_version="nar",    pecan_version="regulated"),
  centrality_cor_results_ex_reg |> mutate(nar_version="nar_ex", pecan_version="regulated")
) |> select(nar_version, pecan_version, everything())

write.csv(
  append_summary(centrality_cor_combined, rho_cen_cols, group_col=c("nar_version","pecan_version")),
  "output/mainstudy/tables/3c_centrality_spearman_combined.csv", row.names = FALSE
)
write.csv(append_summary(centrality_cor_results,     rho_cen_cols),
          "output/mainstudy/tables/3c_centrality_spearman_unreg.csv", row.names = FALSE)
write.csv(append_summary(centrality_cor_results_reg, rho_cen_cols),
          "output/mainstudy/tables/3c_centrality_spearman_reg.csv",   row.names = FALSE)

write.csv(centrality_wilcox_summary,
          "output/mainstudy/tables/3c_centrality_wilcoxon_summary.csv", row.names = FALSE)

# Shared-node degree vectors per VP — for tracing zero-variance exclusions.
# One row per VP × nar_version × shared node.
centrality_shared_nodes <- bind_rows(
  map_dfr(vps, function(vp) {
    shared_ids <- intersect(as.character(nar_degree_list[[vp]]$id),
                            as.character(pecan_nodes_list[[vp]]$id))
    nar_degree_list[[vp]] |>
      mutate(id = as.character(id)) |>
      filter(id %in% shared_ids) |>
      mutate(vp = vp, nar_version = "nar")
  }),
  map_dfr(vps, function(vp) {
    shared_ids <- intersect(as.character(nar_degree_ex_list[[vp]]$id),
                            as.character(pecan_nodes_list[[vp]]$id))
    nar_degree_ex_list[[vp]] |>
      mutate(id = as.character(id)) |>
      filter(id %in% shared_ids) |>
      mutate(vp = vp, nar_version = "nar_ex")
  })
) |>
  select(vp, nar_version, id, in_deg, out_deg)

write.csv(centrality_shared_nodes,
          "output/mainstudy/tables/3c_centrality_shared_node_vectors.csv",
          row.names = FALSE)

# --- 3c-CMP. nar vs. nar_ex centrality comparison ---
# For each Spearman rho column and each PECAN version, compare per-VP rhos
# between nar (full, explicit + implicit) and nar_ex (explicit only).
# Primary comparisons: rho_in_unwt / rho_out_unwt.
# rho_in_wt, rho_out_wt, rho_sevwod retained as sensitivity checks.

centrality_cmp <- bind_rows(
  map_dfr(rho_cen_cols, ~ compare_rho_fn(
    centrality_cor_results[[.x]],     centrality_cor_results_ex[[.x]],
    paste0(.x, " — PECAN unreg: nar vs. nar_ex")
  )),
  map_dfr(rho_cen_cols, ~ compare_rho_fn(
    centrality_cor_results_reg[[.x]], centrality_cor_results_ex_reg[[.x]],
    paste0(.x, " — PECAN reg: nar vs. nar_ex")
  ))
)

cat("\n--- 3b-CMP: nar vs. nar_ex centrality rho comparison ---\n")
print(centrality_cmp)

write.csv(centrality_cmp,
          "output/mainstudy/tables/3c_centrality_nar_vs_ex_comparison.csv",
          row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 3c-TOP. Top-node agreement: does nar identify the same most-central
#             node as PECAN? (indegree and outdegree separately)
# ---------------------------------------------------------------------------
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
# ---------------------------------------------------------------------------

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

# Main function: per-VP top-node agreement for one nar variant vs. one PECAN version
top_node_agree_fn <- function(nar_nodes_lst, nar_edges_lst,
                               pecan_edges_lst, pecan_nodes_lst = pecan_nodes_list) {
  map_dfr(vps, function(vp) {
    nar_deg   <- top_node_fn(nar_edges_lst[[vp]],  nar_nodes_lst[[vp]])
    pecan_deg <- top_node_fn(pecan_edges_lst[[vp]], pecan_nodes_lst[[vp]])

    # All nodes with the maximum degree (ties retained)
    top_nar_in    <- nar_deg   |> filter(in_deg  == max(in_deg))  |> pull(id)
    top_nar_out   <- nar_deg   |> filter(out_deg == max(out_deg)) |> pull(id)
    top_pecan_in  <- pecan_deg |> filter(in_deg  == max(in_deg))  |> pull(id)
    top_pecan_out <- pecan_deg |> filter(out_deg == max(out_deg)) |> pull(id)

    nar_ids <- as.character(nar_nodes_lst[[vp]]$id)

    tibble(
      vp = vp,
      # --- Indegree ---
      top_nar_in_label    = get_lbl(top_nar_in,   nar_nodes_lst[[vp]]),
      top_nar_in_deg      = max(nar_deg$in_deg),
      top_pecan_in_label  = get_lbl(top_pecan_in,  pecan_nodes_lst[[vp]]),
      top_pecan_in_deg    = max(pecan_deg$in_deg),
      pecan_top_in_shared = any(top_pecan_in %in% nar_ids),
      in_agree            = length(intersect(top_nar_in, top_pecan_in)) > 0,
      # --- Outdegree ---
      top_nar_out_label   = get_lbl(top_nar_out,  nar_nodes_lst[[vp]]),
      top_nar_out_deg     = max(nar_deg$out_deg),
      top_pecan_out_label = get_lbl(top_pecan_out, pecan_nodes_lst[[vp]]),
      top_pecan_out_deg   = max(pecan_deg$out_deg),
      pecan_top_out_shared = any(top_pecan_out %in% nar_ids),
      out_agree            = length(intersect(top_nar_out, top_pecan_out)) > 0
    )
  })
}

# Compute for nar and nar_ex vs. unregulated and regulated PECAN
top_agree_nar_unreg    <- top_node_agree_fn(nar_nodes_list,    nar_edges_list,    pecan_edges_list)
top_agree_nar_ex_unreg <- top_node_agree_fn(nar_ex_nodes_list, nar_ex_edges_list, pecan_edges_list)
top_agree_nar_reg      <- top_node_agree_fn(nar_nodes_list,    nar_edges_list,    pecan_edges_reg_list, pecan_nodes_reg_list)
top_agree_nar_ex_reg   <- top_node_agree_fn(nar_ex_nodes_list, nar_ex_edges_list, pecan_edges_reg_list, pecan_nodes_reg_list)

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
      "  %sdegree: agree = %d/%d (%.0f%%) | PECAN top node in nar: %d/%d | agree | shared: %d/%d (%.0f%%)\n",
      measure,
      n_agree, n,
      100 * n_agree / n,
      n_shared, n,
      n_agree_given_shared, n_shared,
      100 * n_agree_given_shared / max(n_shared, 1)
    ))
  }
}

summarise_top_agree(top_agree_nar_unreg,    "nar vs. PECAN unregulated")
summarise_top_agree(top_agree_nar_ex_unreg, "nar_ex vs. PECAN unregulated")
summarise_top_agree(top_agree_nar_reg,      "nar vs. PECAN regulated")
summarise_top_agree(top_agree_nar_ex_reg,   "nar_ex vs. PECAN regulated")

cat("\nPer-participant detail (nar vs. PECAN unregulated):\n")
print(as.data.frame(top_agree_nar_unreg))

cat("\nPer-participant detail (nar_ex vs. PECAN unregulated):\n")
print(as.data.frame(top_agree_nar_ex_unreg))

cat("\nPer-participant detail (nar vs. PECAN regulated):\n")
print(as.data.frame(top_agree_nar_reg))

cat("\nPer-participant detail (nar_ex vs. PECAN regulated):\n")
print(as.data.frame(top_agree_nar_ex_reg))

# Save per-participant tables
write.csv(top_agree_nar_unreg,
          "output/mainstudy/tables/3c_top_node_agreement_nar_unreg.csv",    row.names = FALSE)
write.csv(top_agree_nar_ex_unreg,
          "output/mainstudy/tables/3c_top_node_agreement_nar_ex_unreg.csv", row.names = FALSE)
write.csv(top_agree_nar_reg,
          "output/mainstudy/tables/3c_top_node_agreement_nar_reg.csv",      row.names = FALSE)
write.csv(top_agree_nar_ex_reg,
          "output/mainstudy/tables/3c_top_node_agreement_nar_ex_reg.csv",   row.names = FALSE)

# --- Overview summary table across all four conditions ---
# One row per comparison; columns for in- and out-degree agreement metrics.
#
# Columns:
#   n_in_agree          — participants where nar and PECAN agree on the
#                         highest in-degree node (same node ranked top in both)
#   pct_in_agree        — as percentage of N
#   n_in_shared         — participants where PECAN's top in-degree node
#                         exists anywhere in the nar network (prerequisite for agreement)
#   pct_in_shared       — as percentage of N
#   n_in_agree_given_shared — agreement rate among participants where the
#                             PECAN top node was at least present in nar
#   pct_in_agree_given_shared — as percentage of n_in_shared
#   (same columns repeated for out-degree)

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

top_agree_overview <- make_top_agree_overview(
  list(top_agree_nar_unreg, top_agree_nar_ex_unreg,
       top_agree_nar_reg,   top_agree_nar_ex_reg),
  c("nar vs. PECAN (unregulated)", "nar_ex vs. PECAN (unregulated)",
    "nar vs. PECAN (regulated)",   "nar_ex vs. PECAN (regulated)")
)

cat("\n--- Top-node agreement overview (all four comparisons) ---\n")
print(top_agree_overview)

write.csv(top_agree_overview,
          "output/mainstudy/tables/3c_top_node_agreement_overview.csv",
          row.names = FALSE)

cat("\nSaved top-node agreement tables.\n")

# ---------------------------------------------------------------------------
# --- 3d Node value concordance: nar importance vs. PECAN severity -----
# ---------------------------------------------------------------------------
# nar nodes carry a value ∈ {1 = not most important, 2 = most important}.
# PECAN nodes carry value ∈ {1–10} (participant-rated severity).
# Per-participant Spearman rho over shared nodes (same ID in both networks).
# Run for nar (full) and nar_ex vs. PECAN (unregulated).
# PECAN regulation does not alter node values, so a regulated variant is not
# applicable here (node values are unchanged by edge regulation).

node_val_fn <- function(nar_nodes_lst, label) {
  map_dfr(vps, function(vp) {
    nn <- nar_nodes_lst[[vp]]   |> mutate(id = as.character(id))
    pn <- pecan_nodes_list[[vp]] |> mutate(id = as.character(id))
    shared_ids <- intersect(nn$id, pn$id)
    if (length(shared_ids) < 3)
      return(tibble(vp = vp, nar_version = label,
                    rho = NA_real_, p = NA_real_, n_shared = length(shared_ids)))
    nar_v   <- nn |> filter(id %in% shared_ids) |> arrange(id) |> pull(value)
    pecan_v <- pn |> filter(id %in% shared_ids) |> arrange(id) |> pull(value)
    if (sd(nar_v, na.rm = TRUE) == 0 || sd(pecan_v, na.rm = TRUE) == 0)
      return(tibble(vp = vp, nar_version = label,
                    rho = NA_real_, p = NA_real_, n_shared = length(shared_ids)))
    ct <- cor.test(nar_v, pecan_v, method = "spearman", exact = FALSE)
    tibble(vp = vp, nar_version = label,
           rho = round(ct$estimate, 3), p = round(ct$p.value, 4),
           n_shared = length(shared_ids))
  })
}

node_val_nar    <- node_val_fn(nar_nodes_list,    "nar")
node_val_nar_ex <- node_val_fn(nar_ex_nodes_list, "nar_ex")

print_node_val <- function(res, label) {
  valid_rho <- res$rho[!is.na(res$rho)]
  cat(sprintf("\n%s — Node value concordance (nar importance vs. PECAN severity):\n", label))
  cat(sprintf("  n valid (≥ 3 shared nodes, variation present): %d / %d\n",
              length(valid_rho), nrow(res)))
  if (length(valid_rho) == 0) { cat("  No valid correlations.\n"); return(invisible(NULL)) }
  cat(sprintf("  Mean rho: %.3f | SD: %.3f | Range: %.3f – %.3f\n",
              mean(valid_rho), sd(valid_rho), min(valid_rho), max(valid_rho)))
  if (length(valid_rho) >= 3) {
    wt <- wilcox.test(valid_rho, mu = 0, alternative = "greater", exact = FALSE)
    cat(sprintf("  One-sample Wilcoxon (H1: rho > 0): V = %.0f, p = %.4f\n",
                wt$statistic, wt$p.value))
  }
}

print_node_val(node_val_nar,    "nar")
print_node_val(node_val_nar_ex, "nar_ex")

node_val_combined <- bind_rows(node_val_nar, node_val_nar_ex) |>
  select(nar_version, everything())

write.csv(
  append_summary(node_val_combined, "rho", group_col = "nar_version"),
  "output/mainstudy/tables/3d_node_value_concordance_per_participant.csv",
  row.names = FALSE
)

# --- 3d-CMP. nar vs. nar_ex node value concordance comparison ---

node_val_cmp <- compare_rho_fn(
  node_val_nar$rho, node_val_nar_ex$rho,
  "nar vs. nar_ex — node importance vs. PECAN severity"
)

cat("\n--- 3b-VAL-CMP: nar vs. nar_ex node value rho comparison ---\n")
print(node_val_cmp)

write.csv(node_val_cmp,
          "output/mainstudy/tables/3d_node_value_nar_vs_ex_comparison.csv",
          row.names = FALSE)

# ---------------------------------------------------------------------------
# --- END 3d -------------------------------------------------------------
# ---------------------------------------------------------------------------

# --- 3d-DICH. Dichotomised node value concordance ----------------------
# ---------------------------------------------------------------------------
# nar node values are binary (1 = not most important, 2 = most important).
# To match this scale, PECAN severity (1–10) is dichotomised per participant
# via a median split: nodes with value > median(VP) → 2, otherwise → 1.
# The same per-participant and group-level Spearman analyses are then repeated.

# Per-participant median split of PECAN node severity → {1, 2}
dich_pecan_nodes <- function(pn) {
  med <- median(pn$value, na.rm = TRUE)
  pn |> mutate(value_dich = if_else(value > med, 2L, 1L))
}

# Per-participant Spearman: nar importance vs. dichotomised PECAN severity
node_val_dich_fn <- function(nar_nodes_lst, label) {
  map_dfr(vps, function(vp) {
    nn <- nar_nodes_lst[[vp]]    |> mutate(id = as.character(id))
    pn <- dich_pecan_nodes(pecan_nodes_list[[vp]] |> mutate(id = as.character(id)))
    shared_ids <- intersect(nn$id, pn$id)
    if (length(shared_ids) < 3)
      return(tibble(vp = vp, nar_version = label,
                    rho = NA_real_, p = NA_real_, n_shared = length(shared_ids)))
    nar_v   <- nn  |> filter(id %in% shared_ids) |> arrange(id) |> pull(value)
    pecan_v <- pn  |> filter(id %in% shared_ids) |> arrange(id) |> pull(value_dich)
    if (sd(nar_v, na.rm = TRUE) == 0 || sd(pecan_v, na.rm = TRUE) == 0)
      return(tibble(vp = vp, nar_version = label,
                    rho = NA_real_, p = NA_real_, n_shared = length(shared_ids)))
    ct <- cor.test(nar_v, pecan_v, method = "spearman", exact = FALSE)
    tibble(vp = vp, nar_version = label,
           rho = round(ct$estimate, 3), p = round(ct$p.value, 4),
           n_shared = length(shared_ids))
  })
}

# Save dichotomised PECAN node tables per participant for inspection
dich_pecan_nodes_list <- setNames(
  lapply(vps, function(vp) {
    dich_pecan_nodes(pecan_nodes_list[[vp]] |> mutate(id = as.character(id)))
  }),
  vps
)

write.csv(
  map_dfr(vps, ~ dich_pecan_nodes_list[[.x]] |> mutate(vp = .x)) |> select(vp, everything()),
  "output/mainstudy/tables/3d_pecan_nodes_dichotomised.csv",
  row.names = FALSE
)

node_val_dich_nar    <- node_val_dich_fn(nar_nodes_list,    "nar")
node_val_dich_nar_ex <- node_val_dich_fn(nar_ex_nodes_list, "nar_ex")

print_node_val(node_val_dich_nar,    "nar [PECAN dichotomised]")
print_node_val(node_val_dich_nar_ex, "nar_ex [PECAN dichotomised]")

node_val_dich_combined <- bind_rows(
  node_val_dich_nar,
  node_val_dich_nar_ex
) |> select(nar_version, everything())

write.csv(
  append_summary(node_val_dich_combined, "rho", group_col = "nar_version"),
  "output/mainstudy/tables/3d_node_value_dich_concordance_per_participant.csv",
  row.names = FALSE
)

# Shared-node importance vs. severity vectors per VP — for tracing zero-variance exclusions.
# One row per VP × nar_version × shared node.
node_val_shared_nodes <- bind_rows(
  map_dfr(vps, function(vp) {
    nn <- nar_nodes_list[[vp]] |> mutate(id = as.character(id))
    pn <- dich_pecan_nodes_list[[vp]] |> mutate(id = as.character(id))
    inner_join(
      nn |> select(id, nar_importance = value),
      pn |> select(id, pecan_severity = value, pecan_severity_dich = value_dich),
      by = "id"
    ) |> mutate(vp = vp, nar_version = "nar")
  }),
  map_dfr(vps, function(vp) {
    nn <- nar_ex_nodes_list[[vp]] |> mutate(id = as.character(id))
    pn <- dich_pecan_nodes_list[[vp]] |> mutate(id = as.character(id))
    inner_join(
      nn |> select(id, nar_importance = value),
      pn |> select(id, pecan_severity = value, pecan_severity_dich = value_dich),
      by = "id"
    ) |> mutate(vp = vp, nar_version = "nar_ex")
  })
) |>
  select(vp, nar_version, id, nar_importance, pecan_severity, pecan_severity_dich)

write.csv(node_val_shared_nodes,
          "output/mainstudy/tables/3d_node_value_shared_node_vectors.csv",
          row.names = FALSE)

# --- 3d-DICH-CMP. nar vs. nar_ex dichotomised node value concordance comparison ---

node_val_dich_cmp <- compare_rho_fn(
  node_val_dich_nar$rho, node_val_dich_nar_ex$rho,
  "nar vs. nar_ex — node importance vs. PECAN severity (dichotomised)"
)

cat("\n--- 3d-DICH-CMP: nar vs. nar_ex node value rho comparison (dichotomised) ---\n")
print(node_val_dich_cmp)

write.csv(node_val_dich_cmp,
          "output/mainstudy/tables/3d_node_value_dich_nar_vs_ex_comparison.csv",
          row.names = FALSE)

# ---------------------------------------------------------------------------
# --- END 3d-DICH --------------------------------------------------------
# ---------------------------------------------------------------------------

# ===========================================================================
# SECTION 4: DESCRIPTIVE & AGGREGATED DATA
# ===========================================================================
#relevant for the master thesis were the following output tables:
#4a_ext_anova_network_metrics #RM-ANOVA: tests descriptive network metrics for significant differences within a metric between network types
#4a_ext_posthoc_network_metrics #Bonferroni corrected paired t-tests: tests descriptive network metrics for significant differences within a metric between network types
#group_network_nar
#group_network_nar_ex
#group_network_PECAN_reg
#group_network_PECAN_unreg

cat("\n=== SECTION 4: DESCRIPTIVE & AGGREGATED DATA ===\n")
# Note: all_nar_nodes/edges and all_pecan_nodes/edges were pooled in Section 0.



# --- 4a-EXT. Network metrics compared across 4 network types (RM-ANOVA + post-hoc) ---
# Metrics: n nodes, n edges, n feedback loops (≤4), density
# Networks: nar (full), nar_ex, PECAN (unreg), PECAN (reg)
# Design: one-way within-subjects (repeated measures) ANOVA, 4 levels, n = 50 VPs.
# Sphericity: Greenhouse-Geisser correction applied automatically by afex.
# Post-hoc: all pairwise Bonferroni-corrected t-tests via emmeans::pairs().

# Build long-format metric table (one row per VP × network)
metric_data_long <- map_dfr(vps, function(vp) {
  ne    <- distinct(nar_edges_list[[vp]],    from, to)
  ne_ex <- distinct(nar_ex_edges_list[[vp]], from, to)
  pe    <- pecan_edges_list[[vp]]
  pe_r  <- pecan_edges_reg_list[[vp]]

  nn    <- nrow(nar_nodes_list[[vp]])
  nn_ex <- nrow(nar_ex_nodes_list[[vp]])
  np    <- nrow(pecan_nodes_list[[vp]])
  np_r  <- nrow(pecan_nodes_reg_list[[vp]])

  dens <- function(n_e, n_n) if (n_n > 1) n_e / (n_n * (n_n - 1)) else NA_real_

  loops_nar    <- count_cycles_le4(build_adj(ne,    nar_nodes_list[[vp]]$id))$total_le4
  loops_nar_ex <- count_cycles_le4(build_adj(ne_ex, nar_ex_nodes_list[[vp]]$id))$total_le4
  loops_pecan  <- count_cycles_le4(build_adj(pe,    pecan_nodes_list[[vp]]$id))$total_le4
  loops_pecan_r <- count_cycles_le4(build_adj(pe_r, pecan_nodes_reg_list[[vp]]$id))$total_le4

  bind_rows(
    tibble(vp = vp, network = "nar",         n_nodes = nn,    n_edges = nrow(ne),    n_loops = loops_nar,     density = dens(nrow(ne),    nn)),
    tibble(vp = vp, network = "nar_ex",      n_nodes = nn_ex, n_edges = nrow(ne_ex), n_loops = loops_nar_ex,  density = dens(nrow(ne_ex), nn_ex)),
    tibble(vp = vp, network = "PECAN_unreg", n_nodes = np,    n_edges = nrow(pe),    n_loops = loops_pecan,   density = dens(nrow(pe),    np)),
    tibble(vp = vp, network = "PECAN_reg",   n_nodes = np_r,  n_edges = nrow(pe_r),  n_loops = loops_pecan_r, density = dens(nrow(pe_r),  np_r))
  )
}) |>
  mutate(
    network = factor(network, levels = c("nar", "nar_ex", "PECAN_unreg", "PECAN_reg")),
    vp      = factor(vp)
  )

write.csv(metric_data_long,
          "output/mainstudy/tables/4a_ext_network_metrics_long.csv", row.names = FALSE)

# Helper: RM-ANOVA + Bonferroni post-hoc for one metric
rm_anova_ph_fn <- function(metric_col, metric_label) {
  cat(sprintf("\n=== %s ===\n", metric_label))
  fit <- afex::aov_ez(
    id          = "vp",
    dv          = metric_col,
    data        = metric_data_long,
    within      = "network",
    anova_table = list(es = "pes", correction = "GG")
  )
  cat("RM-ANOVA (GG-corrected):\n"); print(fit)
  emm <- emmeans::emmeans(fit, ~ network)
  ph  <- pairs(emm, adjust = "bonferroni")
  cat("Post-hoc pairwise (Bonferroni):\n"); print(ph)
  list(anova = fit, posthoc = ph)
}

cat("\n--- 4a-EXT: RM-ANOVA — n nodes, n edges, n loops, density across network types ---\n")
anova_n_nodes <- rm_anova_ph_fn("n_nodes", "N Nodes")
anova_n_edges <- rm_anova_ph_fn("n_edges", "N Edges")
anova_n_loops <- rm_anova_ph_fn("n_loops", "N Feedback Loops (≤4)")
anova_density <- rm_anova_ph_fn("density", "Density")

# Save ANOVA tables and post-hoc results to CSV
anova_results_list <- list(
  n_nodes = anova_n_nodes, n_edges = anova_n_edges,
  n_loops = anova_n_loops, density = anova_density
)
anova_tbl <- map_dfr(names(anova_results_list), \(m)
  as.data.frame(anova_results_list[[m]]$anova$anova_table) |>
    mutate(metric = m, .before = 1)
)
posthoc_tbl <- map_dfr(names(anova_results_list), \(m)
  as.data.frame(summary(anova_results_list[[m]]$posthoc)) |>
    mutate(metric = m, .before = 1)
)
write.csv(anova_tbl,   "output/mainstudy/tables/4a_ext_anova_network_metrics.csv",   row.names = FALSE)
write.csv(posthoc_tbl, "output/mainstudy/tables/4a_ext_posthoc_network_metrics.csv", row.names = FALSE)
cat("[4a-EXT complete]\n")


# ---------------------------------------------------------------------------
# --- 4b. SIF: Symptom Importance Factor ---
# Formula: SIF = proportional_frequency × (mean_severity / max_scale)
# PECAN max scale = 10 | nar has no severity, so SIF_nar = proportional frequency only
# ---------------------------------------------------------------------------

sif_pecan <- all_pecan_nodes |>
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
    prop_freq = round(n_reported / n_vps, 3),
    SIF       = round(prop_freq * (mean_severity / 10), 3)
  ) |>
  arrange(desc(SIF))

sif_nar_fn <- function(all_nodes, suffix) {
  all_nodes |>
    group_by(id, label) |>
    summarise(n_reported = n(), .groups = "drop") |>
    mutate(prop_freq = round(n_reported / n_vps, 3),
           !!paste0("SIF_", suffix) := prop_freq) |>
    arrange(desc(prop_freq))
}

sif_nar    <- sif_nar_fn(all_nar_nodes,    "nar")
sif_nar_ex <- sif_nar_fn(all_nar_ex_nodes, "nar_ex")

cat("Top 10 symptoms by SIF (PECAN):\n");  print(head(sif_pecan, 10))
cat("\nTop 10 symptoms by SIF_nar:\n");     print(head(sif_nar,    10))
cat("\nTop 10 symptoms by SIF_nar_ex:\n");  print(head(sif_nar_ex, 10))

write.csv(sif_pecan,  "output/mainstudy/tables/4b_SIF_PECAN.csv",  row.names = FALSE)
write.csv(sif_nar,    "output/mainstudy/tables/4b_SIF_nar.csv",    row.names = FALSE)
write.csv(sif_nar_ex, "output/mainstudy/tables/4b_SIF_nar_ex.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 4c. EIF: Edge Impact Factor ---
# Formula: EIF = (n_VPs_with_edge / n_VPs_both_nodes_present) × (mean_width / max_scale)
# PECAN max_scale = 5 | nar all widths = 1 so EIF_nar = proportional frequency
# ---------------------------------------------------------------------------

# Pre-compute PECAN co-occurrence (nar co-occurrence built later in 5b with ex/im)
pecan_cooccur        <- build_cooccurrence(pecan_nodes_list, vps)
pecan_cooccur_counts <- pecan_cooccur |> count(from, to, name = "n_both_present")

# PECAN edge frequencies
pecan_edge_freq <- all_pecan_edges |>
  group_by(from, to) |>
  summarise(
    n_edge_reported = n(),
    mean_width      = round(mean(width, na.rm = TRUE), 2),
    sd_width        = round(sd(width,   na.rm = TRUE), 2),
    min_width       = min(width, na.rm = TRUE),
    max_width       = max(width, na.rm = TRUE),
    .groups = "drop"
  )

eif_pecan <- pecan_edge_freq |>
  left_join(pecan_cooccur_counts, by = c("from", "to")) |>
  mutate(
    prop_freq = if_else(
      !is.na(n_both_present) & n_both_present > 0,
      round(n_edge_reported / n_both_present, 3),
      NA_real_
    ),
    EIF = round(prop_freq * (mean_width / 5), 3)
  ) |>
  left_join(all_pecan_nodes |> select(id, from_label = label) |> distinct(),
            by = c("from" = "id")) |>
  left_join(all_pecan_nodes |> select(id, to_label = label) |> distinct(),
            by = c("to" = "id")) |>
  select(from_label, to_label, from, to, n_edge_reported, n_both_present,
         prop_freq, mean_width, sd_width, min_width, max_width, EIF) |>
  arrange(desc(EIF))

# nar / nar_ex edge frequencies and EIF
eif_nar_fn <- function(all_edges, all_nodes_sub, cooccur_counts, suffix) {
  edge_freq <- all_edges |>
    distinct(from, to, vp) |>
    group_by(from, to) |>
    summarise(n_edge_reported = n(), .groups = "drop")
  edge_freq |>
    left_join(cooccur_counts, by = c("from", "to")) |>
    mutate(
      prop_freq = if_else(!is.na(n_both_present) & n_both_present > 0,
                          round(n_edge_reported / n_both_present, 3), NA_real_),
      !!paste0("EIF_", suffix) := prop_freq
    ) |>
    left_join(all_nodes_sub |> select(id, from_label=label) |> distinct(), by=c("from"="id")) |>
    left_join(all_nodes_sub |> select(id, to_label=label)   |> distinct(), by=c("to"="id")) |>
    select(from_label, to_label, from, to, n_edge_reported, n_both_present,
           prop_freq, starts_with("EIF_")) |>
    arrange(desc(prop_freq))
}

nar_cooccur    <- build_cooccurrence(nar_nodes_list,    vps)
nar_ex_cooccur <- build_cooccurrence(nar_ex_nodes_list, vps)

nar_cooccur_counts    <- nar_cooccur    |> count(from, to, name = "n_both_present")
nar_ex_cooccur_counts <- nar_ex_cooccur |> count(from, to, name = "n_both_present")

nar_edge_freq <- all_nar_edges |> distinct(from, to, vp) |>
  group_by(from, to) |> summarise(n_edge_reported = n(), .groups = "drop")

eif_nar    <- eif_nar_fn(all_nar_edges,    all_nar_nodes,    nar_cooccur_counts,    "nar")
eif_nar_ex <- eif_nar_fn(all_nar_ex_edges, all_nar_ex_nodes, nar_ex_cooccur_counts, "nar_ex")

cat("Top 10 edges by EIF (PECAN):\n");  print(head(eif_pecan,  10))
cat("\nTop 10 edges by EIF_nar:\n");    print(head(eif_nar,    10))
cat("\nTop 10 edges by EIF_nar_ex:\n"); print(head(eif_nar_ex, 10))

write.csv(eif_pecan,  "output/mainstudy/tables/4c_EIF_PECAN_unreg.csv", row.names = FALSE)
write.csv(eif_nar,    "output/mainstudy/tables/4c_EIF_nar.csv",         row.names = FALSE)
write.csv(eif_nar_ex, "output/mainstudy/tables/4c_EIF_nar_ex.csv",      row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 4c (REGULATED). EIF on regulated PECAN networks ---
# Co-occurrence counts are unchanged (node-based); only edge frequencies differ.
# ---------------------------------------------------------------------------

pecan_edge_freq_reg <- all_pecan_edges_reg |>
  group_by(from, to) |>
  summarise(
    n_edge_reported = n(),
    mean_width      = round(mean(width, na.rm = TRUE), 2),
    sd_width        = round(sd(width,   na.rm = TRUE), 2),
    min_width       = min(width, na.rm = TRUE),
    max_width       = max(width, na.rm = TRUE),
    .groups = "drop"
  )

eif_pecan_reg <- pecan_edge_freq_reg |>
  left_join(pecan_cooccur_counts, by = c("from", "to")) |>
  mutate(
    prop_freq = if_else(
      !is.na(n_both_present) & n_both_present > 0,
      round(n_edge_reported / n_both_present, 3),
      NA_real_
    ),
    EIF = round(prop_freq * (mean_width / 5), 3)
  ) |>
  left_join(all_pecan_nodes |> select(id, from_label = label) |> distinct(),
            by = c("from" = "id")) |>
  left_join(all_pecan_nodes |> select(id, to_label = label) |> distinct(),
            by = c("to" = "id")) |>
  select(from_label, to_label, from, to, n_edge_reported, n_both_present,
         prop_freq, mean_width, sd_width, min_width, max_width, EIF) |>
  arrange(desc(EIF))

cat("\n[REGULATED] Top 10 edges by EIF (PECAN):\n")
print(head(eif_pecan_reg, 10))

write.csv(eif_pecan_reg, "output/mainstudy/tables/4c_EIF_PECAN_reg.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 4d. Feedback Loops (matrix algebra, cycles up to length 4) ---
# Count directed simple cycles of length 2, 3, and 4 using trace of adjacency
# matrix powers. This is O(n³) and bounded regardless of network size.
## For binary directed adjacency matrix A (no self-loops):
#   n_2cycles = trace(A²) / 2
#   n_3cycles = trace(A³) / 3
#   n_4cycles = (trace(A⁴) − 2·Σᵢ(A²[i,i])² + trace(A²)) / 4
#
# The n_4cycles correction subtracts non-simple length-4 closed walks:
#   Case v₁=v₃ (wings node revisited): contributes Σᵢ (A²[i,i])² to trace(A⁴)
#   Case v₂=v₀ (start revisited mid-walk): contributes Σᵢ (A²[i,i])² to trace(A⁴)
#   Overlap (both cases simultaneously, e.g. i→j→i→j→i): contributes trace(A²)
# Inclusion-exclusion: 2·Σ(A²diag²) − trace(A²), subtracted from trace(A⁴).
# ---------------------------------------------------------------------------

# PECAN loop counts
loop_pecan <- map_dfr(vps, function(vp) {
  adj <- build_adj(pecan_edges_list[[vp]], pecan_nodes_list[[vp]]$id)
  cyc <- count_cycles_le4(adj)
  tibble(vp = vp,
         pecan_n2  = cyc$n2,  pecan_n3  = cyc$n3,
         pecan_n4  = cyc$n4,  pecan_le4 = cyc$total_le4)
})

# nar / nar_ex loop counts (edges binary; deduplicate within VP)
loop_nar <- map_dfr(vps, function(vp) {
  adj <- build_adj(nar_edges_list[[vp]], nar_nodes_list[[vp]]$id)
  cyc <- count_cycles_le4(adj)
  tibble(vp=vp, nar_n2=cyc$n2, nar_n3=cyc$n3, nar_n4=cyc$n4, nar_le4=cyc$total_le4)
})
loop_nar_ex <- map_dfr(vps, function(vp) {
  adj <- build_adj(nar_ex_edges_list[[vp]], nar_ex_nodes_list[[vp]]$id)
  cyc <- count_cycles_le4(adj)
  tibble(vp=vp, nar_ex_n2=cyc$n2, nar_ex_n3=cyc$n3, nar_ex_n4=cyc$n4, nar_ex_le4=cyc$total_le4)
})

# Regulated PECAN loop counts
loop_pecan_reg <- map_dfr(vps, function(vp) {
  adj <- build_adj(pecan_edges_reg_list[[vp]], pecan_nodes_list[[vp]]$id)
  cyc <- count_cycles_le4(adj)
  tibble(vp = vp,
         pecan_reg_n2  = cyc$n2,  pecan_reg_n3  = cyc$n3,
         pecan_reg_n4  = cyc$n4,  pecan_reg_le4 = cyc$total_le4)
})

loop_results <- left_join(loop_pecan, loop_nar, by = "vp") |>
  left_join(loop_pecan_reg, by = "vp") |>
  left_join(loop_nar_ex, by = "vp")

cat("\nFeedback loops per participant (PECAN unreg | nar | PECAN reg | nar_ex, cycles ≤ 4):\n")
print(loop_results)

cat("\n--- Summary ---\n")
cat("PECAN unregulated  mean loops (≤4):", round(mean(loop_results$pecan_le4),     1),
    "| SD:", round(sd(loop_results$pecan_le4),     1), "\n")
cat("PECAN regulated    mean loops (≤4):", round(mean(loop_results$pecan_reg_le4), 1),
    "| SD:", round(sd(loop_results$pecan_reg_le4), 1), "\n")
cat("nar                mean loops (≤4):", round(mean(loop_results$nar_le4),       1),
    "| SD:", round(sd(loop_results$nar_le4),       1), "\n")
cat("nar_ex             mean loops (≤4):", round(mean(loop_results$nar_ex_le4),    1),
    "| SD:", round(sd(loop_results$nar_ex_le4),    1), "\n")

# Paired test: nar vs. PECAN unregulated
loop_diff_unreg <- loop_results$nar_le4 - loop_results$pecan_le4
sw_loops_unreg  <- shapiro.test(loop_diff_unreg)
cat("\nShapiro-Wilk (nar − PECAN unreg diff): W =", round(sw_loops_unreg$statistic, 3),
    ", p =", round(sw_loops_unreg$p.value, 4), "\n")
cat("\nPaired test (nar vs. PECAN unregulated, n_loops_le4):\n")
if (sw_loops_unreg$p.value > 0.05) {
  print(t.test(loop_results$nar_le4, loop_results$pecan_le4, paired = TRUE))
} else {
  print(wilcox.test(loop_results$nar_le4, loop_results$pecan_le4, paired = TRUE))
}

# Paired test: nar vs. PECAN regulated
loop_diff_reg <- loop_results$nar_le4 - loop_results$pecan_reg_le4
sw_loops_reg  <- shapiro.test(loop_diff_reg)
cat("\nShapiro-Wilk (nar − PECAN reg diff): W =", round(sw_loops_reg$statistic, 3),
    ", p =", round(sw_loops_reg$p.value, 4), "\n")
cat("\n[REGULATED] Paired test (nar vs. PECAN regulated, n_loops_le4):\n")
if (sw_loops_reg$p.value > 0.05) {
  print(t.test(loop_results$nar_le4, loop_results$pecan_reg_le4, paired = TRUE))
} else {
  print(wilcox.test(loop_results$nar_le4, loop_results$pecan_reg_le4, paired = TRUE))
}

write.csv(loop_results, "output/mainstudy/tables/4d_feedback_loops.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 4d-ii. Common loops across participants ---
# Enumerate all simple directed cycles up to length 4 by node label,
# then count how many participants share each cycle.
#
# Canonical form: lexicographically smallest cyclic rotation of the label
# sequence (direction preserved; only the starting node is shifted).
# ---------------------------------------------------------------------------


# --- PECAN common loops ---
pecan_cycle_list <- setNames(
  lapply(vps, function(vp)
    enumerate_cycles_labeled(pecan_edges_list[[vp]], pecan_nodes_list[[vp]])),
  vps
)

pecan_common_loops <- pecan_cycle_list |>
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
cat("Total unique cycles found:", nrow(pecan_common_loops), "\n")
cat("Cycles shared by ≥ 2 participants:\n")
print(pecan_common_loops |> filter(n_participants >= 2))

# --- nar common loops ---
nar_cycle_list <- setNames(
  lapply(vps, function(vp)
    enumerate_cycles_labeled(nar_edges_list[[vp]], nar_nodes_list[[vp]])),
  vps
)

nar_common_loops <- nar_cycle_list |>
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
cat("Total unique cycles found:", nrow(nar_common_loops), "\n")
cat("Cycles shared by ≥ 2 participants:\n")
print(nar_common_loops |> filter(n_participants >= 2))

write.csv(pecan_common_loops, "output/mainstudy/tables/4d_common_loops_PECAN.csv", row.names = FALSE)
write.csv(nar_common_loops,   "output/mainstudy/tables/4d_common_loops_nar.csv",   row.names = FALSE)

# nar_ex common loops
nar_ex_cycle_list <- setNames(
  lapply(vps, function(vp)
    enumerate_cycles_labeled(nar_ex_edges_list[[vp]], nar_ex_nodes_list[[vp]])),
  vps)
nar_ex_common_loops <- nar_ex_cycle_list |>
  imap_dfr(~ tibble(vp=.y, cycle=.x)) |>
  filter(!is.na(cycle)) |>
  group_by(cycle) |>
  summarise(n_participants=n(), participants=paste(sort(vp), collapse=", "), .groups="drop") |>
  mutate(cycle_length=str_count(cycle, "\u2192")+1L) |>
  arrange(desc(n_participants), cycle_length, cycle)
cat("\nCommon feedback loops — nar_ex (cycles ≤ 4):\n")
cat("Total unique:", nrow(nar_ex_common_loops), "| Shared by ≥ 2:\n")
print(nar_ex_common_loops |> filter(n_participants >= 2))
write.csv(nar_ex_common_loops, "output/mainstudy/tables/4d_common_loops_nar_ex.csv", row.names=FALSE)



# ---------------------------------------------------------------------------
# --- 4e. Aggregated Group Networks ---
# Inclusion threshold: ≥ 10% of participants (computed dynamically from n_vps)
# ---------------------------------------------------------------------------

threshold <- ceiling(0.10 * n_vps)

# PECAN group network (UNREGULATED)
pecan_node_include <- sif_pecan |> filter(n_reported >= threshold)
pecan_edge_include <- eif_pecan |>
  filter(n_edge_reported >= threshold,
         from %in% pecan_node_include$id,
         to   %in% pecan_node_include$id)

cat("[UNREGULATED] PECAN group network: ", nrow(pecan_node_include), "nodes,",
    nrow(pecan_edge_include), "edges\n")

# matching the style used in loop_visualisation_mainstudy.R.
plot_grp_network_vis <- function(nodes_df, edges_df, title, color, outfile) {
  net <- visNetwork(nodes_df, edges_df) |>
    visNodes(shape = "dot", size = 20,
             color = list(background = "#2E86AB", border = "#2E86AB",
                          highlight   = list(background = "#2E86AB", border = "#2E86AB")),
             font  = list(color = "black", size = 16,
                          strokeWidth = 4, strokeColor = "white")) |>
    visEdges(arrows = "to",
             color  = list(color = color, highlight = color)) |>
    visIgraphLayout(layout = "layout_with_fr") |>
    visOptions(manipulation = FALSE)
  tmp <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(net, tmp, selfcontained = TRUE)
  webshot2::webshot(tmp, file = outfile, vwidth = 1600, vheight = 1200, delay = 2)
  cat("Saved:", outfile, "\n")
}

if (nrow(pecan_edge_include) > 0 && nrow(pecan_node_include) > 0) {
  grp_ids  <- unique(c(pecan_edge_include$from, pecan_edge_include$to))
  grp_info <- pecan_node_include |> filter(id %in% grp_ids) |> distinct(id, .keep_all = TRUE)
  plot_grp_network_vis(
    nodes_df = data.frame(id = grp_info$id, label = grp_info$label),
    edges_df  = data.frame(from = pecan_edge_include$from, to = pecan_edge_include$to),
    title    = paste0("PECAN Group Network UNREGULATED (≥", threshold, "/", n_vps, " VPs)"),
    color    = "#2E86AB",
    outfile  = "output/mainstudy/plots/group_network_PECAN_unreg.png"
  )
}

# PECAN group network (REGULATED)
pecan_node_include_reg <- sif_pecan |> filter(n_reported >= threshold)
pecan_edge_include_reg <- eif_pecan_reg |>
  filter(n_edge_reported >= threshold,
         from %in% pecan_node_include_reg$id,
         to   %in% pecan_node_include_reg$id)

cat("[REGULATED] PECAN group network: ", nrow(pecan_node_include_reg), "nodes,",
    nrow(pecan_edge_include_reg), "edges\n")

if (nrow(pecan_edge_include_reg) > 0 && nrow(pecan_node_include_reg) > 0) {
  grp_ids_reg  <- unique(c(pecan_edge_include_reg$from, pecan_edge_include_reg$to))
  grp_info_reg <- pecan_node_include_reg |> filter(id %in% grp_ids_reg) |> distinct(id, .keep_all = TRUE)
  plot_grp_network_vis(
    nodes_df = data.frame(id = grp_info_reg$id, label = grp_info_reg$label),
    edges_df  = data.frame(from = pecan_edge_include_reg$from, to = pecan_edge_include_reg$to),
    title    = paste0("PECAN Group Network REGULATED (≥", threshold, "/", n_vps, " VPs)"),
    color    = "#2E86AB",
    outfile  = "output/mainstudy/plots/group_network_PECAN_reg.png"
  )
}

# Helper to build and plot a nar-variant group network
plot_nar_grp_network <- function(sif_tbl, eif_tbl, label, color, file_suffix) {
  n_inc <- sif_tbl |> filter(n_reported      >= threshold)
  e_inc <- eif_tbl |> filter(n_edge_reported >= threshold,
                               from %in% n_inc$id, to %in% n_inc$id)
  cat(sprintf("%s group network: %d nodes, %d edges\n", label, nrow(n_inc), nrow(e_inc)))
  if (nrow(e_inc) == 0 || nrow(n_inc) == 0) return(invisible(NULL))
  gids  <- unique(c(e_inc$from, e_inc$to))
  ginfo <- n_inc |> filter(id %in% gids)
  plot_grp_network_vis(
    nodes_df = data.frame(id = ginfo$id, label = ginfo$label),
    edges_df  = data.frame(from = e_inc$from, to = e_inc$to),
    title    = paste0(label, " Group Network (≥", threshold, "/", n_vps, " VPs)"),
    color    = color,
    outfile  = sprintf("output/mainstudy/plots/group_network_%s.png", file_suffix)
  )
}

# nar group network
nar_node_include <- sif_nar |> filter(n_reported >= threshold)
nar_edge_include <- eif_nar |>
  filter(n_edge_reported >= threshold,
         from %in% nar_node_include$id, to %in% nar_node_include$id)
plot_nar_grp_network(sif_nar, eif_nar, "nar", "#2E86AB", "nar")

# nar_ex group network
plot_nar_grp_network(sif_nar_ex, eif_nar_ex, "nar_ex", "#2E86AB", "nar_ex")


# ---------------------------------------------------------------------------
# --- 4e-ii. Aggregated group network summary metrics ---
# For each group network (nar, nar_ex, PECAN unregulated, PECAN regulated)
# report: n_nodes and n_edges surviving the ≥ 10% threshold, network density,
# and directed feedback loop counts (cycles of length 2–4).
# ---------------------------------------------------------------------------

grp_net_metrics <- map_dfr(
  list(
    list(label = "nar",        sif = sif_nar,    eif = eif_nar),
    list(label = "nar_ex",     sif = sif_nar_ex, eif = eif_nar_ex),
    list(label = "PECAN_unreg",sif = sif_pecan,  eif = eif_pecan),
    list(label = "PECAN_reg",  sif = sif_pecan,  eif = eif_pecan_reg)
  ),
  function(nw) {
    n_inc   <- nw$sif |> filter(n_reported      >= threshold)
    e_inc   <- nw$eif |> filter(n_edge_reported >= threshold,
                                 from %in% n_inc$id,
                                 to   %in% n_inc$id)
    n_nodes <- nrow(n_inc)
    n_edges <- nrow(e_inc)
    density <- if (n_nodes > 1) round(n_edges / (n_nodes * (n_nodes - 1)), 4) else NA_real_
    adj     <- build_adj(e_inc, as.character(n_inc$id))
    cyc     <- count_cycles_le4(adj)
    tibble(
      network     = nw$label,
      n_nodes     = n_nodes,
      n_edges     = n_edges,
      density     = density,
      n_2cycles   = cyc$n2,
      n_3cycles   = cyc$n3,
      n_4cycles   = cyc$n4,
      n_loops_le4 = cyc$total_le4
    )
  }
)

cat("\nAggregated group network metrics (threshold ≥", threshold, "/", n_vps, "VPs):\n")
print(grp_net_metrics)

write.csv(grp_net_metrics,
          "output/mainstudy/tables/4e_aggregated_network_metrics.csv",
          row.names = FALSE)

# ---------------------------------------------------------------------------
# --- 4f. Symptom frequency table (combined nar + nar_ex + PECAN) ---
# ---------------------------------------------------------------------------

symptom_table <- full_join(
  sif_nar    |> select(id, label,    n_nar    = n_reported, pct_nar    = prop_freq),
  sif_nar_ex |> select(id, ex_label = label,  n_nar_ex = n_reported, pct_nar_ex = prop_freq),
  by = "id"
) |>
  full_join(sif_pecan |> select(id, pecan_label=label, n_pecan=n_reported,
                                 pct_pecan=prop_freq, mean_severity, sd_severity,
                                 min_severity, max_severity, SIF),
            by = "id") |>
  mutate(
    label    = coalesce(label, ex_label, pecan_label),
    n_nar    = replace_na(n_nar,    0L),
    n_nar_ex = replace_na(n_nar_ex, 0L),
    n_pecan  = replace_na(n_pecan,  0L)
  ) |>
  select(-ex_label, -pecan_label) |>
  arrange(desc(n_pecan + n_nar))
# ---------------------------------------------------------------------------
# --- 4g. nar vs. PECAN edge comparison ---
# How often was each directed edge reported in nar vs. PECAN?
# ---------------------------------------------------------------------------

# Label lookup combining both networks
label_lookup <- bind_rows(
  all_nar_nodes   |> select(id, label),
  all_pecan_nodes |> select(id, label)
) |>
  distinct(id, .keep_all = TRUE)

edge_comparison <- full_join(
  nar_edge_freq   |> select(from, to, n_nar = n_edge_reported),
  pecan_edge_freq |> select(from, to, n_pecan = n_edge_reported),
  by = c("from", "to")
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
  left_join(label_lookup |> rename(from_label = label), by = c("from" = "id")) |>
  left_join(label_lookup |> rename(to_label   = label), by = c("to"   = "id")) |>
  select(from_label, to_label, from, to, n_nar, n_pecan, diff, more_in) |>
  arrange(desc(abs(diff)))

cat("\nEdge comparison — most divergent (nar vs. PECAN):\n")
print(head(edge_comparison, 20))
cat("\nSummary of edge direction:\n")
print(count(edge_comparison, more_in))

write.csv(edge_comparison, "output/mainstudy/tables/4g_edge_comparison_nar_vs_PECAN_unreg.csv",
          row.names = FALSE)

# --- 4h (REGULATED). nar vs. PECAN regulated edge comparison ---

edge_comparison_reg <- full_join(
  nar_edge_freq        |> select(from, to, n_nar   = n_edge_reported),
  pecan_edge_freq_reg  |> select(from, to, n_pecan = n_edge_reported),
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
  left_join(label_lookup |> rename(from_label = label), by = c("from" = "id")) |>
  left_join(label_lookup |> rename(to_label   = label), by = c("to"   = "id")) |>
  select(from_label, to_label, from, to, n_nar, n_pecan, diff, more_in) |>
  arrange(desc(abs(diff)))

cat("\n[REGULATED] Edge comparison — most divergent (nar vs. PECAN regulated):\n")
print(head(edge_comparison_reg, 20))
cat("\n[REGULATED] Summary of edge direction:\n")
print(count(edge_comparison_reg, more_in))

write.csv(edge_comparison_reg, "output/mainstudy/tables/4g_edge_comparison_nar_vs_PECAN_reg.csv",
          row.names = FALSE)

# ===========================================================================
# SECTION 5: SUBJECTIVE REPRESENTATIVENESS (3-WAY RM-ANOVA)
# ===========================================================================
#relevant for the master thesis were the following output tables and plots:
#5b_representativeness_descriptives_reduced #descriptives of subjective representativeness and feasibility ratings
#5b_representativeness_anova_reduced #ANOVA comparing the subjective representativeness of the networks (nar_full, nar_ex, PECAN reg)
#representativeness_profile_plot_reduced.png

# Participants rated representativeness of their networks after viewing them.
# Three conditions (within-subjects):
#   NW02 = PECAN  |  NW06 = nar (full: ex + im)  |  NW07 = nar_ex
#
# Reversed items (higher score = LESS representative → flip before summing):
#   NW02_14, NW02_15, NW06_14, NW06_15, NW07_14, NW07_15  (6 - x, scale 1–5)
# t2_scored (with all reversals applied) was built in Section 2e.

# --- 5.1 T2 data ---
# t2 was loaded in Section 2a-T2; t2_scored (all 6 items reversed) from Section 2e.
cat("\nT2 data ready:", nrow(t2), "participants\n")

# --- 5.2 Item lists ---

nw02_items <- names(t2_scored)[grepl("^NW02_", names(t2_scored))]
nw06_items <- names(t2_scored)[grepl("^NW06_", names(t2_scored))]
nw07_items <- names(t2_scored)[grepl("^NW07_", names(t2_scored))]

cat("NW02 items (PECAN,  n =", length(nw02_items), "):", paste(nw02_items, collapse = ", "), "\n")
cat("NW06 items (nar,    n =", length(nw06_items), "):", paste(nw06_items, collapse = ", "), "\n")
cat("NW07 items (nar_ex, n =", length(nw07_items), "):", paste(nw07_items, collapse = ", "), "\n")
cat("Reversed: NW02_14, NW02_15, NW06_14, NW06_15, NW07_14, NW07_15 (applied in Section 2e)\n")

# --- 5.3 Sum scores ---

rep_data <- t2_scored |>
  select(VP, all_of(c(nw02_items, nw06_items, nw07_items))) |>
  mutate(
    PECAN  = rowSums(across(all_of(nw02_items)), na.rm = FALSE),
    nar    = rowSums(across(all_of(nw06_items)), na.rm = FALSE),
    nar_ex = rowSums(across(all_of(nw07_items)), na.rm = FALSE)
  ) |>
  select(VP, PECAN, nar, nar_ex)

cat("\nSum scores per participant:\n")
print(rep_data)
cat("\nMissing values: PECAN =", sum(is.na(rep_data$PECAN)),
    "| nar =",    sum(is.na(rep_data$nar)),
    "| nar_ex =", sum(is.na(rep_data$nar_ex)), "\n")

# Complete cases: all three conditions present
rep_complete <- rep_data |> filter(!is.na(PECAN) & !is.na(nar) & !is.na(nar_ex))
cat("Complete cases for RM-ANOVA: n =", nrow(rep_complete), "\n")

# --- 5.4 Descriptives ---

rep_long <- rep_complete |>
  pivot_longer(cols = c(PECAN, nar, nar_ex),
               names_to  = "model",
               values_to = "representativeness") |>
  mutate(VP    = factor(VP),
         model = factor(model, levels = c("PECAN", "nar_ex", "nar")))

max_score <- length(nw02_items) * 5

rep_desc <- rep_long |>
  group_by(model) |>
  summarise(
    M    = round(mean(representativeness), 2),
    SD   = round(sd(representativeness),   2),
    Min  = min(representativeness),
    Max  = max(representativeness),
    .groups = "drop"
  ) |>
  mutate(
    M_pct   = round(M   / max_score * 100, 1),
    SD_pct  = round(SD  / max_score * 100, 1),
    Min_pct = round(Min / max_score * 100, 1),
    Max_pct = round(Max / max_score * 100, 1)
  )

cat("\nDescriptives (sum scores, max =", max_score, "| percentages of max):\n")
print(rep_desc)

# --- 5.5 Assumption checks ---
# Normality: Shapiro-Wilk on each condition's distribution.
# (With small n — see nrow(rep_complete) — power for Shapiro-Wilk is low; inspect alongside descriptives.)

cat("\nShapiro-Wilk normality checks:\n")
for (cond in c("PECAN", "nar_ex", "nar")) {
  x  <- rep_complete[[cond]]
  sw <- shapiro.test(x)
  cat(sprintf("  %-8s W = %.3f, p = %.4f\n", cond, sw$statistic, sw$p.value))
}

# --- 5.6 One-way within-subjects RM-ANOVA (3 conditions) ---
# afex::aov_ez applies Greenhouse-Geisser correction when sphericity is violated
# (Mauchly's test W and ε are reported automatically in the output).

rep_anova <- aov_ez(
  id          = "VP",
  dv          = "representativeness",
  data        = rep_long,
  within      = "model",
  anova_table = list(es = "pes")
)

cat("\nRM-ANOVA — subjective representativeness (PECAN vs. nar_ex vs. nar):\n")
print(rep_anova)

# --- 5.7 Post-hoc pairwise comparisons (Bonferroni-corrected) ---
# Only run if the omnibus F is significant (α = .05); reported regardless for completeness.

emm_rep <- emmeans(rep_anova, ~ model)
rep_posthoc <- pairs(emm_rep, adjust = "bonferroni")

cat("\nPost-hoc pairwise comparisons (Bonferroni-corrected):\n")
print(rep_posthoc)

# Effect sizes (Cohen's d) for each pair
rep_posthoc_tbl <- as.data.frame(summary(rep_posthoc)) |>
  rename(pair = contrast, t_ratio = t.ratio, p_bonf = p.value)

# --- 5.7b Non-parametric sensitivity: Friedman test ---
# Used as robustness check if RM-ANOVA assumptions are questionable.

rep_mat <- rep_complete |>
  select(PECAN, nar_ex, nar) |>
  as.matrix()

friedman_rep <- friedman.test(rep_mat)
cat("\nFriedman test (non-parametric sensitivity):\n")
print(friedman_rep)

# Post-hoc Wilcoxon pairwise if Friedman significant
wilcox_pairs <- list(
  PECAN_vs_nar_ex = wilcox.test(rep_complete$PECAN,  rep_complete$nar_ex, paired=TRUE, exact=FALSE),
  PECAN_vs_nar    = wilcox.test(rep_complete$PECAN,  rep_complete$nar,    paired=TRUE, exact=FALSE),
  nar_ex_vs_nar   = wilcox.test(rep_complete$nar_ex, rep_complete$nar,    paired=TRUE, exact=FALSE)
)
wilcox_tbl <- map_dfr(names(wilcox_pairs), \(nm) tibble(
  pair = nm,
  V    = wilcox_pairs[[nm]]$statistic,
  p_uncorrected = round(wilcox_pairs[[nm]]$p.value, 4),
  p_bonferroni  = round(min(wilcox_pairs[[nm]]$p.value * 3, 1), 4)
))
cat("\nWilcoxon signed-rank pairwise (Bonferroni-corrected p × 3):\n")
print(wilcox_tbl)

# --- 5.7c Profile plot ---

ggplot(rep_long, aes(x = model, y = representativeness, group = VP)) +
  geom_line(colour = "grey70", linewidth = 0.4) +
  geom_point(colour = "grey50", size = 2) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               colour = "#2E86AB", linewidth = 1.3, linetype = "dashed") +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               shape = 18, size = 5, colour = "#2E86AB") +
  scale_x_discrete(labels = c("PECAN", "Nar-Ex", "Nar")) +
  scale_y_continuous(limits = c(0, max_score), breaks = seq(0, max_score, 10)) +
  labs(
    x = "Network model",
    y = "Representativeness and feasibility (sum score)"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title   = element_text(face = "bold", size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 13),
        axis.text    = element_text(size = 14))

ggsave("output/mainstudy/plots/representativeness_profile_plot.png",
       width = 7, height = 5, dpi = 150)

# --- 5.8 Save outputs ---

write.csv(rep_data,
          "output/mainstudy/tables/5a_representativeness_sumscores.csv",    row.names = FALSE)
write.csv(rep_desc,
          "output/mainstudy/tables/5a_representativeness_descriptives.csv", row.names = FALSE)
write.csv(as.data.frame(rep_anova$anova_table),
          "output/mainstudy/tables/5a_representativeness_anova.csv",        row.names = TRUE)
write.csv(rep_posthoc_tbl,
          "output/mainstudy/tables/5a_representativeness_posthoc.csv",      row.names = FALSE)
write.csv(wilcox_tbl,
          "output/mainstudy/tables/5a_representativeness_wilcox_pairs.csv", row.names = FALSE)

cat("Saved representativeness tables.\n")

# ===========================================================================
# SECTION 5b: SUBJECTIVE REPRESENTATIVENESS — NETWORK-SPECIFIC ITEMS ONLY
# ===========================================================================
# Items NW02_05–NW02_09 (and the equivalent NW06_05–NW06_09, NW07_05–NW07_09)
# ask about networks in general rather than about the specific network shown.
# They are therefore excluded here to obtain a sum score that reflects only
# network-specific representativeness judgements.
#
# Everything else is identical to Section 5 (same reversal, same ANOVA, same
# plots/tables). Variables carry the suffix _red ("reduced item set").
# Reversed items: same as Section 5; NW*_14 and NW*_15 are reversed in
# t2_scored already (applied in Section 2d).

# --- 5b.2 Reduced item lists (exclude items 05–09 from each scale) ---

general_items <- paste0("0", 5:9)   # "05","06","07","08","09"

nw02_items_red <- nw02_items[!sub("^NW02_", "", nw02_items) %in% general_items]
nw06_items_red <- nw06_items[!sub("^NW06_", "", nw06_items) %in% general_items]
nw07_items_red <- nw07_items[!sub("^NW07_", "", nw07_items) %in% general_items]

cat("\n[6b] Network-specific items only (NW02_05–NW02_09 excluded):\n")
cat("NW02 items (PECAN,  n =", length(nw02_items_red), "):", paste(nw02_items_red, collapse = ", "), "\n")
cat("NW06 items (nar,    n =", length(nw06_items_red), "):", paste(nw06_items_red, collapse = ", "), "\n")
cat("NW07 items (nar_ex, n =", length(nw07_items_red), "):", paste(nw07_items_red, collapse = ", "), "\n")
cat("Excluded (general network items):",
    paste(paste0("NW02_", general_items), collapse = ", "), "and equivalents for NW06/NW07\n")

# --- 5b.3 Sum scores (reduced) ---

rep_data_red <- t2_scored |>
  select(VP, all_of(c(nw02_items_red, nw06_items_red, nw07_items_red))) |>
  mutate(
    PECAN  = rowSums(across(all_of(nw02_items_red)), na.rm = FALSE),
    nar    = rowSums(across(all_of(nw06_items_red)), na.rm = FALSE),
    nar_ex = rowSums(across(all_of(nw07_items_red)), na.rm = FALSE)
  ) |>
  select(VP, PECAN, nar, nar_ex)

cat("\n[6b] Reduced sum scores per participant:\n")
print(rep_data_red)
cat("\nMissing values: PECAN =", sum(is.na(rep_data_red$PECAN)),
    "| nar =",    sum(is.na(rep_data_red$nar)),
    "| nar_ex =", sum(is.na(rep_data_red$nar_ex)), "\n")

rep_complete_red <- rep_data_red |> filter(!is.na(PECAN) & !is.na(nar) & !is.na(nar_ex))
cat("Complete cases for RM-ANOVA: n =", nrow(rep_complete_red), "\n")

# --- 5b.4 Descriptives ---

rep_long_red <- rep_complete_red |>
  pivot_longer(cols = c(PECAN, nar, nar_ex),
               names_to  = "model",
               values_to = "representativeness") |>
  mutate(VP    = factor(VP),
         model = factor(model, levels = c("PECAN", "nar_ex", "nar")))

max_score_red <- length(nw02_items_red) * 5

rep_desc_red <- rep_long_red |>
  group_by(model) |>
  summarise(
    M    = round(mean(representativeness), 2),
    SD   = round(sd(representativeness),   2),
    Min  = min(representativeness),
    Max  = max(representativeness),
    .groups = "drop"
  ) |>
  mutate(
    M_pct   = round(M   / max_score_red * 100, 1),
    SD_pct  = round(SD  / max_score_red * 100, 1),
    Min_pct = round(Min / max_score_red * 100, 1),
    Max_pct = round(Max / max_score_red * 100, 1)
  )

cat("\n[5b] Descriptives (reduced sum scores, max =", max_score_red, "| percentages of max):\n")
print(rep_desc_red)

# --- 5b.5 Assumption checks ---

cat("\n[5b] Shapiro-Wilk normality checks (reduced scores):\n")
for (cond in c("PECAN", "nar_ex", "nar")) {
  x  <- rep_complete_red[[cond]]
  sw <- shapiro.test(x)
  cat(sprintf("  %-8s W = %.3f, p = %.4f\n", cond, sw$statistic, sw$p.value))
}

# --- 5b.6 One-way within-subjects RM-ANOVA (reduced scores) ---

rep_anova_red <- aov_ez(
  id          = "VP",
  dv          = "representativeness",
  data        = rep_long_red,
  within      = "model",
  anova_table = list(es = "pes")
)

cat("\n[6b] RM-ANOVA — representativeness, network-specific items only:\n")
print(rep_anova_red)

# --- 5b.7 Post-hoc pairwise comparisons (Bonferroni-corrected) ---

emm_rep_red     <- emmeans(rep_anova_red, ~ model)
rep_posthoc_red <- pairs(emm_rep_red, adjust = "bonferroni")

cat("\n[6b] Post-hoc pairwise comparisons (Bonferroni-corrected):\n")
print(rep_posthoc_red)

rep_posthoc_red_tbl <- as.data.frame(summary(rep_posthoc_red)) |>
  rename(pair = contrast, t_ratio = t.ratio, p_bonf = p.value)

# --- 5b.7b Non-parametric sensitivity: Friedman test ---

rep_mat_red <- rep_complete_red |>
  select(PECAN, nar_ex, nar) |>
  as.matrix()

friedman_rep_red <- friedman.test(rep_mat_red)
cat("\n[5b] Friedman test (non-parametric sensitivity):\n")
print(friedman_rep_red)

wilcox_pairs_red <- list(
  PECAN_vs_nar_ex = wilcox.test(rep_complete_red$PECAN,  rep_complete_red$nar_ex, paired=TRUE, exact=FALSE),
  PECAN_vs_nar    = wilcox.test(rep_complete_red$PECAN,  rep_complete_red$nar,    paired=TRUE, exact=FALSE),
  nar_ex_vs_nar   = wilcox.test(rep_complete_red$nar_ex, rep_complete_red$nar,    paired=TRUE, exact=FALSE)
)
wilcox_tbl_red <- map_dfr(names(wilcox_pairs_red), \(nm) tibble(
  pair = nm,
  V    = wilcox_pairs_red[[nm]]$statistic,
  p_uncorrected = round(wilcox_pairs_red[[nm]]$p.value, 4),
  p_bonferroni  = round(min(wilcox_pairs_red[[nm]]$p.value * 3, 1), 4)
))
cat("\n[5b] Wilcoxon signed-rank pairwise (Bonferroni-corrected p × 3):\n")
print(wilcox_tbl_red)

# --- 5b.7c Profile plot (reduced scores) ---

ggplot(rep_long_red, aes(x = model, y = representativeness, group = VP)) +
  geom_line(colour = "grey70", linewidth = 0.4) +
  geom_point(colour = "grey50", size = 2) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               colour = "#2E86AB", linewidth = 1.3, linetype = "dashed") +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               shape = 18, size = 5, colour = "#2E86AB") +
  scale_x_discrete(labels = c("PECAN", "Nar-Ex", "Nar")) +
  scale_y_continuous(limits = c(0, max_score_red), breaks = seq(0, max_score_red, 5)) +
  labs(
    x = "Network model",
    y = "Representativeness and feasibility (sum score)"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title   = element_text(face = "bold", size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 13),
        axis.text    = element_text(size = 14))

ggsave("output/mainstudy/plots/representativeness_profile_plot_reduced.png",
       width = 7, height = 5, dpi = 150)

# --- 5b.8 Save outputs ---

write.csv(rep_data_red,
          "output/mainstudy/tables/5b_representativeness_sumscores_reduced.csv",    row.names = FALSE)
write.csv(rep_desc_red,
          "output/mainstudy/tables/5b_representativeness_descriptives_reduced.csv", row.names = FALSE)
write.csv(as.data.frame(rep_anova_red$anova_table),
          "output/mainstudy/tables/5b_representativeness_anova_reduced.csv",        row.names = TRUE)
write.csv(rep_posthoc_red_tbl,
          "output/mainstudy/tables/5b_representativeness_posthoc_reduced.csv",      row.names = FALSE)
write.csv(wilcox_tbl_red,
          "output/mainstudy/tables/5b_representativeness_wilcox_pairs_reduced.csv", row.names = FALSE)

cat("Saved reduced-item representativeness tables (Section 5b).\n")

# --- Save all key results to a single RData file ---

save(
  # ── Section 1: demographics ──────────────────────────────────────────────
  t1, age_stats, gender_table, education_table,
  # ── Section 2: DASS / reliability ────────────────────────────────────────
  dass_data, alpha_result,
  t2, dass_t2_data,
  alpha_nw_tbl,           # Cronbach α for NW02 / NW06 / NW07
  jaccard_results, jaccard_summary,
  jacc_nar_unreg, jacc_nar_reg, jacc_ex_unreg, jacc_ex_reg,
  pecan_nodes_reg_list,
  jacc_nar_unreg, jacc_nar_reg, jacc_ex_unreg, jacc_ex_reg,
  jacc_pair_tests,
  # ── Section 3a: edge concordance — nar / nar_ex × unreg / reg ────────────
  edge_cor_results,    edge_cor_results_ex,
  edge_cor_results_reg, edge_cor_results_ex_reg,
  edge_cor_combined,
  # ── Section 3b: per-participant centrality ────────────────────────────────
  centrality_cor_results,    centrality_cor_results_ex,
  centrality_cor_results_reg, centrality_cor_results_ex_reg,
  centrality_cor_combined,
  # ── Section 3a-CMP: edge concordance ex vs im ─────────────────────────────
  edge_concord_cmp,
  # ── Section 3b-CMP: centrality ex vs im ──────────────────────────────────
  centrality_cmp,
  # ── Section 3b-VAL: node value concordance ──────────────────────────────
  node_val_nar, node_val_nar_ex, node_val_combined,
  # ── Section 3b-VAL-CMP / 3d-DICH-CMP: node value comparisons ────────────
  node_val_cmp, node_val_dich_cmp,
  # ── Section 3b-VAL-DICH: dichotomised PECAN severity ─────────────────────
  node_val_dich_nar, node_val_dich_nar_ex, node_val_dich_combined,
  # ── Section 3c: node / edge counts ───────────────────────────────────────

  # ── Section 4: descriptives ───────────────────────────────────────────────
  sif_pecan, sif_nar, sif_nar_ex,
  eif_pecan, eif_pecan_reg, eif_nar, eif_nar_ex,
  loop_results,
  pecan_common_loops, nar_common_loops, nar_ex_common_loops,
  grp_net_metrics,
  symptom_table,
  edge_comparison, edge_comparison_reg,
  # ── Section 5: representativeness ────────────────────────────────────────
  rep_data, rep_desc,
  rep_anova, rep_posthoc_tbl, friedman_rep, wilcox_tbl,
  # ── Section 5b: reduced-item representativeness (NW*_05–09 excluded) ─────
  rep_data_red, rep_desc_red,
  rep_anova_red, rep_posthoc_red_tbl, friedman_rep_red, wilcox_tbl_red,
  file = "output/mainstudy/N2N_results.RData"
)

# ===========================================================================
# SECTION 6: REPORTING OVERVIEW TABLES
# ===========================================================================
#relevant for the master thesis were the following output tables:
#6_overview_01_sample_descriptives
#6_overview_02_network_descriptives

# Compact summary tables for thesis reporting (mean, SD, min, max).

cat("\n=== SECTION 6: REPORTING OVERVIEW TABLES ===\n")

sumstat <- function(x, label, digits = 2) {
  x <- x[!is.na(x)]
  tibble(label = label, n = length(x),
         M = round(mean(x), digits), SD = round(sd(x), digits),
         Min = round(min(x), digits), Max = round(max(x), digits))
}
# ── Table 1: Sample descriptives (DASS-21) ───────────────────────────────────
tbl_sample <- bind_rows(
  sumstat(dass_data$DASS_total,       "DASS-21 total score — T1 (0–63)"),
  sumstat(dass_t2_data$DASS_T2_total, "DASS-21 total score — T2 (0–63)")
)
cat("\n[Table 1] Sample descriptives:\n"); print(tbl_sample)
write.csv(tbl_sample, "output/mainstudy/tables/6_overview_01_sample_descriptives.csv",
          row.names = FALSE)

# ── Table 2: Network descriptives ────────────────────────────────────────────
# Node/edge counts, density, severity, and mean edge width per network version.

pecan_n_nodes        <- map_int(pecan_nodes_list,     nrow)
pecan_n_nodes_reg    <- map_int(pecan_nodes_reg_list, nrow)
nar_n_nodes          <- map_int(nar_nodes_list,        nrow)
nar_ex_n_nodes       <- map_int(nar_ex_nodes_list,     nrow)
pecan_sev_per_vp     <- map_dbl(pecan_nodes_list,     ~ mean(.x$value, na.rm = TRUE))
pecan_sev_per_vp_reg <- map_dbl(pecan_nodes_reg_list, ~ mean(.x$value, na.rm = TRUE))
pecan_mean_width     <- map_dbl(pecan_edges_list,     ~ mean(.x$width, na.rm = TRUE))
pecan_mean_width_reg <- map_dbl(pecan_edges_reg_list, ~ mean(.x$width, na.rm = TRUE))

# Directed-graph density: edges / (n * (n - 1)); returns NA when n < 2.
.dens <- function(n_edges, n_nodes) {
  ifelse(n_nodes < 2, NA_real_, n_edges / (n_nodes * (n_nodes - 1)))
}

pecan_n_edges_v      <- map_int(vps, ~ nrow(pecan_edges_list[[.x]]))
pecan_strength_v     <- map_dbl(vps, ~ sum(pecan_edges_list[[.x]]$width,     na.rm = TRUE))
pecan_density_v      <- .dens(pecan_n_edges_v,     pecan_n_nodes)
pecan_n_edges_reg_v  <- map_int(vps, ~ nrow(pecan_edges_reg_list[[.x]]))
pecan_strength_reg_v <- map_dbl(vps, ~ sum(pecan_edges_reg_list[[.x]]$width, na.rm = TRUE))
pecan_density_reg_v  <- .dens(pecan_n_edges_reg_v, pecan_n_nodes_reg)
nar_n_edges_v        <- map_int(vps, ~ nrow(distinct(nar_edges_list[[.x]],    from, to)))
nar_density_v        <- .dens(nar_n_edges_v,    nar_n_nodes)
nar_ex_n_edges_v     <- map_int(vps, ~ nrow(distinct(nar_ex_edges_list[[.x]], from, to)))
nar_ex_density_v     <- .dens(nar_ex_n_edges_v, nar_ex_n_nodes)

tbl_networks <- bind_rows(
  # PECAN unregulated
  sumstat(pecan_n_nodes,        "PECAN (unreg) — nodes per participant"),
  sumstat(pecan_n_edges_v,      "PECAN (unreg) — edges per participant"),
  sumstat(pecan_strength_v,     "PECAN (unreg) — total edge weight (strength)"),
  sumstat(pecan_mean_width,     "PECAN (unreg) — mean edge width (1–5)"),
  sumstat(pecan_density_v,      "PECAN (unreg) — density"),
  sumstat(pecan_sev_per_vp,     "PECAN (unreg) — mean node severity (1–10)"),
  # PECAN regulated
  sumstat(pecan_n_nodes_reg,    "PECAN (reg)   — nodes per participant"),
  sumstat(pecan_n_edges_reg_v,  "PECAN (reg)   — edges per participant"),
  sumstat(pecan_strength_reg_v, "PECAN (reg)   — total edge weight (strength)"),
  sumstat(pecan_mean_width_reg, "PECAN (reg)   — mean edge width (1–5)"),
  sumstat(pecan_density_reg_v,  "PECAN (reg)   — density"),
  sumstat(pecan_sev_per_vp_reg, "PECAN (reg)   — mean node severity (1–10)"),
  # nar (full: explicit + implicit)
  sumstat(nar_n_nodes,          "nar — nodes per participant"),
  sumstat(nar_n_edges_v,        "nar — edges per participant"),
  sumstat(nar_density_v,        "nar — density"),
  # nar_ex (explicit only)
  sumstat(nar_ex_n_nodes,       "nar_ex — nodes per participant"),
  sumstat(nar_ex_n_edges_v,     "nar_ex — edges per participant"),
  sumstat(nar_ex_density_v,     "nar_ex — density"),
  # Feedback loops (cycles ≤ 4)
  sumstat(loop_results$pecan_le4,    "PECAN (unreg) — feedback loops (≤4) per participant"),
  sumstat(loop_results$pecan_reg_le4,"PECAN (reg)   — feedback loops (≤4) per participant"),
  sumstat(loop_results$nar_le4,      "nar           — feedback loops (≤4) per participant"),
  sumstat(loop_results$nar_ex_le4,   "nar_ex        — feedback loops (≤4) per participant")
)
cat("\n[Table 2] Network descriptives:\n"); print(tbl_networks, n = Inf)
write.csv(tbl_networks, "output/mainstudy/tables/6_overview_02_network_descriptives.csv",
          row.names = FALSE)

# ===========================================================================
# SECTION 7: RAIN PLOTS — SPEARMAN CORRELATION DISTRIBUTIONS
# ===========================================================================
#relevant for the master thesis were the following plots:
#dotplot_N_Nodes
#dotplot_N_Edges
#dotplot_N_Feedback_Loops
#dotplot_Density
#rainplot_node_severity
#rainplot_edge_concordance
#rainplot_centrality_outdegree
#rainplot_centrality_indegree

# One plot per metric: edges, node severity, indegree, outdegree.
# X-axis: comparison condition (nar version × PECAN regulation).
# Fill:   PECAN version (unregulated vs regulated).
# Node severity has no regulation dimension → 3 groups only, faceted by
# severity type (original 1–10 vs per-participant median-split dichotomised).
# ---------------------------------------------------------------------------

library(ggdist)

# Colour palette: unregulated = blue, regulated = orange
pecan_colours <- c("unregulated" = "#909090", "regulated" = "#2E86AB")

# Factor levels for rain plots — internal keys keep (unreg)/(reg) for correct
# ordering and colour mapping; scale_x_discrete in rain_layers() strips them
# for display since colour already encodes the PECAN regulation type.
comparison_levels <- c(
  "Nar\n(unreg)", "Nar\n(reg)",
  "Nar-Ex\n(unreg)", "Nar-Ex\n(reg)"
)

make_comparison <- function(df, nar_col = "nar_version", pecan_col = "pecan_version") {
  nar_label_map <- c("nar" = "Nar", "nar_ex" = "Nar-Ex")
  df |>
    mutate(
      nar_label = nar_label_map[.data[[nar_col]]],
      comparison = factor(
        paste0(nar_label, "\n(", sub("ulated", "", .data[[pecan_col]]), ")"),
        levels = comparison_levels
      ),
      pecan_version = factor(.data[[pecan_col]],
                             levels = c("unregulated", "regulated"))
    ) |>
    select(-nar_label)
}

# Core raincloud layer helper (returns a list of ggplot layers)
rain_layers <- function(fill_var = "pecan_version") {
  list(
    ggdist::stat_halfeye(
      aes(fill = .data[[fill_var]]),
      adjust = 0.6, width = 0.5, .width = 0,      justification = -0.25, point_colour = NA, alpha = 0.75
    ),
    geom_boxplot(
      aes(fill = .data[[fill_var]]),
      width = 0.12, outlier.shape = NA, alpha = 0.6
    ),
    ggdist::stat_dots(
      aes(fill = .data[[fill_var]], colour = .data[[fill_var]]),
      side = "left", justification = 1.12,
      dotsize = 0.6, binwidth = 0.04, alpha = 0.8
    ),
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40"),
    coord_cartesian(ylim = c(-1.05, 1.05)),
    # Strip (unreg)/(reg) from x-axis display; colour legend already encodes it
    scale_x_discrete(labels = c(
      "Nar\n(unreg)"   = "Nar",
      "Nar\n(reg)"     = "Nar",
      "Nar-Ex\n(unreg)" = "Nar-Ex",
      "Nar-Ex\n(reg)"   = "Nar-Ex"
    )),
    labs(y = "Spearman \u03c1", x = NULL,
         fill = "PECAN", colour = "PECAN"),
    theme_minimal(base_size = 14),
    theme(legend.position = "bottom",
          axis.text.x  = element_text(size = 14, lineheight = 0.85),
          axis.title.y = element_text(size = 16),
          plot.title   = element_text(face = "bold", size = 14))
  )
}

# ---------------------------------------------------------------------------

edge_plot_df <- read.csv(
  "output/mainstudy/tables/3b_edge_concordance_spearman_combined.csv"
) |>
  filter(!is.na(rho)) |>
  make_comparison()

p_edges <- ggplot(edge_plot_df, aes(x = comparison, y = rho)) +
  rain_layers() +
  scale_fill_manual(values = pecan_colours) +
  scale_colour_manual(values = pecan_colours) +
  labs(y = "Spearman \u03c1", x = NULL)

ggsave("output/mainstudy/plots/rainplot_edge_concordance.png",
       p_edges, width = 8, height = 5, dpi = 300)
cat("\n[Section 7a] Edge concordance rain plot saved.\n")

# ---------------------------------------------------------------------------
# 7b. Node severity concordance (original + dichotomised, faceted)
# ---------------------------------------------------------------------------

# Original (1–10 PECAN severity)
nv_orig <- node_val_combined |>
  filter(!is.na(rho)) |>
  mutate(
    nar_version   = factor(nar_version, levels = c("nar", "nar_ex")),
    severity_type = "Original (1\u201310)"
  )

# Dichotomised (per-participant median split)
nv_dich <- node_val_dich_combined |>
  filter(!is.na(rho)) |>
  mutate(
    nar_version   = factor(nar_version, levels = c("nar", "nar_ex")),
    severity_type = "Dichotomised (median split)"
  )

nv_plot_df <- bind_rows(nv_orig, nv_dich) |>
  mutate(severity_type = factor(severity_type,
                                levels = c("Original (1\u201310)",
                                           "Dichotomised (median split)")))

p_node_sev <- ggplot(
  nv_plot_df |>
    filter(severity_type == "Dichotomised (median split)") |>
    mutate(nar_version = recode(nar_version, "nar" = "Nar", "nar_ex" = "Nar-Ex")),
  aes(x = nar_version, y = rho)
) +
  ggdist::stat_halfeye(
    adjust = 0.6, width = 0.5, .width = 0,
    justification = -0.25, point_colour = NA, alpha = 0.75,
    fill = "#2E86AB"
  ) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6,
               fill = "#2E86AB") +
  ggdist::stat_dots(
    side = "left", justification = 1.12,
    dotsize = 0.6, binwidth = 0.04, alpha = 0.8,
    fill = "#909090", colour = "#909090"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  coord_cartesian(ylim = c(-1.05, 1.05)) +
  labs(y = "Spearman \u03c1", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title   = element_text(face = "bold", size = 14),
        axis.text.x  = element_text(size = 14),
        axis.title.y = element_text(size = 16))

ggsave("output/mainstudy/plots/rainplot_node_severity.png",
       p_node_sev, width = 8, height = 5, dpi = 300)
cat("[Section 7b] Node severity concordance rain plot saved.\n")

# ---------------------------------------------------------------------------
# 7c & 7d. Centrality: indegree and outdegree (unweighted, primary measure)
# ---------------------------------------------------------------------------

cen_plot_df <- read.csv(
  "output/mainstudy/tables/3c_centrality_spearman_combined.csv"
) |>
  filter(n_shared >= 3) |>
  make_comparison()

# Indegree (unweighted)
p_indegree <- ggplot(
  cen_plot_df |> filter(!is.na(rho_in_unwt)),
  aes(x = comparison, y = rho_in_unwt)
) +
  rain_layers() +
  scale_fill_manual(values = pecan_colours) +
  scale_colour_manual(values = pecan_colours) +
  labs(y = "Spearman \u03c1", x = NULL)

ggsave("output/mainstudy/plots/rainplot_centrality_indegree.png",
       p_indegree, width = 8, height = 5, dpi = 300)
cat("[Section 7c] Indegree centrality rain plot saved.\n")

# Outdegree (unweighted)
p_outdegree <- ggplot(
  cen_plot_df |> filter(!is.na(rho_out_unwt)),
  aes(x = comparison, y = rho_out_unwt)
) +
  rain_layers() +
  scale_fill_manual(values = pecan_colours) +
  scale_colour_manual(values = pecan_colours) +
  labs(y = "Spearman \u03c1", x = NULL)

ggsave("output/mainstudy/plots/rainplot_centrality_outdegree.png",
       p_outdegree, width = 8, height = 5, dpi = 300)
cat("[Section 7d] Outdegree centrality rain plot saved.\n")

cat("\n=== Section 7a–7d complete. Rain plots saved to output/mainstudy/plots/ ===\n")

# ---------------------------------------------------------------------------
# 7e. Dot plots — one plot per metric, panels = network version
# ---------------------------------------------------------------------------
# Each plot shows the per-participant distribution of one metric across all
# applicable network versions (nar, nar_ex, PECAN unreg, PECAN reg).
# Metrics that only exist in PECAN (Strength, Mean Node Severity, Mean Edge
# Width) show 2 panels; shared metrics show all 4.
# Blue diamond = mean; blue bars = ±1 SD.

# Sort loop_results to match vps order
lr <- loop_results[match(vps, loop_results$vp), ]

network_levels <- c("Nar", "Nar-Ex", "PECAN (unreg)", "PECAN (reg)")
metric_levels  <- c(
  "N Nodes", "N Edges", "N Feedback Loops", "Density",
  "Strength (Weighted)",
  "Mean Node Severity (1-10)", "Mean Edge Width (1-5)"
)

# Build combined long-format data; missing metrics for nar/nar_ex stay NA
# Edge counts, strength, and density use the inline vectors computed in Section 6.
all_dot_long <- bind_rows(
  tibble(
    vp = vps, network = "Nar",
    `N Nodes`          = as.double(nar_n_nodes[vps]),
    `N Edges`          = as.double(nar_n_edges_v),
    `N Feedback Loops` = lr$nar_le4,
    `Density`          = nar_density_v
  ),
  tibble(
    vp = vps, network = "Nar-Ex",
    `N Nodes`          = as.double(nar_ex_n_nodes[vps]),
    `N Edges`          = as.double(nar_ex_n_edges_v),
    `N Feedback Loops` = lr$nar_ex_le4,
    `Density`          = nar_ex_density_v
  ),
  tibble(
    vp = vps, network = "PECAN (unreg)",
    `N Nodes`                   = as.double(pecan_n_nodes[vps]),
    `N Edges`                   = as.double(pecan_n_edges_v),
    `N Feedback Loops`          = lr$pecan_le4,
    `Density`                   = pecan_density_v,
    `Strength (Weighted)`       = pecan_strength_v,
    `Mean Node Severity (1-10)` = pecan_sev_per_vp[vps],
    `Mean Edge Width (1-5)`     = pecan_mean_width[vps]
  ),
  tibble(
    vp = vps, network = "PECAN (reg)",
    `N Nodes`                   = as.double(pecan_n_nodes[vps]),
    `N Edges`                   = as.double(pecan_n_edges_reg_v),
    `N Feedback Loops`          = lr$pecan_reg_le4,
    `Density`                   = pecan_density_reg_v,
    `Strength (Weighted)`       = pecan_strength_reg_v,
    `Mean Node Severity (1-10)` = pecan_sev_per_vp[vps],
    `Mean Edge Width (1-5)`     = pecan_mean_width_reg[vps]
  )
) |>
  pivot_longer(-c(vp, network), names_to = "metric", values_to = "value") |>
  mutate(
    network = factor(network, levels = network_levels),
    metric  = factor(metric,  levels = metric_levels)
  ) |>
  filter(!is.na(value))

# Piecewise linear y-axis transformation.
# The axis is split at break_val: the lower segment [0, break_val] occupies
# break_frac of the total axis height; the upper segment [break_val, max_val]
# occupies the remainder. This keeps the denser lower range readable while
# still showing the full extent of outlying networks.
piecewise_trans <- function(break_val, break_frac, max_val) {
  trans_fn <- function(x) {
    ifelse(x <= break_val,
           break_frac * x / break_val,
           break_frac + (1 - break_frac) * (x - break_val) / (max_val - break_val))
  }
  inv_fn <- function(y) {
    ifelse(y <= break_frac,
           y * break_val / break_frac,
           break_val + (y - break_frac) * (max_val - break_val) / (1 - break_frac))
  }
  scales::trans_new("piecewise_linear", trans_fn, inv_fn)
}

# Three-segment piecewise linear transformation.
# Each segment [0, b1], [b1, b2], [b2, max_val] occupies exactly 1/3 of axis height.
three_seg_trans <- function(b1, b2, max_val) {
  trans_fn <- function(x) {
    ifelse(x <= b1,
      x / b1 / 3,
      ifelse(x <= b2,
        1/3 + (x - b1) / (b2 - b1) / 3,
        2/3 + (x - b2) / (max_val - b2) / 3
      )
    )
  }
  inv_fn <- function(y) {
    ifelse(y <= 1/3,
      y * 3 * b1,
      ifelse(y <= 2/3,
        b1 + (y - 1/3) * 3 * (b2 - b1),
        b2 + (y - 2/3) * 3 * (max_val - b2)
      )
    )
  }
  scales::trans_new("three_seg_linear", trans_fn, inv_fn)
}

# Pre-compute maxima needed for the transformations.
# Limits are rounded up to clean numbers so the single upper-range label is tidy.
# Only one break is placed in the compressed upper segment to prevent overlap.
max_loops_dot <- max(all_dot_long |> filter(metric == "N Feedback Loops") |> pull(value), na.rm = TRUE)
max_edges_dot <- max(all_dot_long |> filter(metric == "N Edges")          |> pull(value), na.rm = TRUE)

max_loops_lim <- ceiling(max_loops_dot / 100) * 100   # e.g. 1938 → 2000
max_edges_lim <- ceiling(max_edges_dot / 10)  * 10    # e.g. 138  →  140

# Breaks chosen so each of the three axis segments has clear reference labels.
# 30→100 gap in transformed space ≈ 0.42 in; all others ≥ 0.9 in.
loops_axis_breaks <- c(0, 10, 20, 30, 100, 250, max_loops_lim)
edges_axis_breaks <- c(0, 10, 25, 50,      max_edges_lim)

# One plot per metric
metric_dot_plot <- function(metric_name) {
  df       <- all_dot_long |> filter(metric == metric_name)
  n_panels <- n_distinct(df$network)

  p <- ggplot(df, aes(x = 0, y = value)) +
    geom_jitter(width = 0.18, height = 0, size = 1.8, alpha = 0.6,
                colour = "grey50") +
    stat_summary(
      fun.data = mean_sdl, fun.args = list(mult = 1),
      geom = "errorbar", width = 0.35,
      colour = "#2E86AB", linewidth = 0.6
    ) +
    stat_summary(
      fun = mean, geom = "point",
      size = 3.5, shape = 18, colour = "#2E86AB"
    ) +
    facet_wrap(~ network, nrow = 1) +
    labs(x = NULL, y = metric_name) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y       = element_text(size = 16),
      strip.text         = element_text(size = 14, face = "bold"),
      plot.title         = element_text(face = "bold", size = 14),
      panel.spacing      = unit(0.8, "lines")
    )

  # Custom piecewise y-axis for metrics dominated by one network type.
  # Single label in the upper (compressed) segment avoids overlap.
  if (metric_name == "N Feedback Loops") {
    # Three equal thirds: 0–30 | 30–250 | 250–max_loops_lim
    p <- p + scale_y_continuous(
      trans  = three_seg_trans(30, 250, max_loops_lim),
      breaks = loops_axis_breaks,
      limits = c(0, max_loops_lim),
      oob    = scales::squish
    )
  } else if (metric_name == "N Edges") {
    # Lower 50% of axis: 0–50; upper 50%: 50–max_edges_lim
    p <- p + scale_y_continuous(
      trans  = piecewise_trans(50, 0.5, max_edges_lim),
      breaks = edges_axis_breaks,
      limits = c(0, max_edges_lim),
      oob    = scales::squish
    )
  }

  p
}

# Save one PNG per metric
walk(levels(all_dot_long$metric), function(m) {
  p     <- metric_dot_plot(m)
  fname <- paste0("dotplot_", gsub("[^a-zA-Z0-9]+", "_", m), ".png")
  df_m  <- all_dot_long |> filter(metric == m)
  w     <- n_distinct(df_m$network) * 2.4 + 0.6
  ggsave(paste0("output/mainstudy/plots/", fname),
         p, width = w, height = 4, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
})

cat("[Section 7e] Per-metric dot plots saved.\n")

cat("\n=== Section 7 complete. Rain plots and dot plots saved to output/mainstudy/plots/ ===\n")

cat("\n=== Analysis complete. All results saved to output/ ===\n")
cat("Tables: output/mainstudy/tables/\n")
cat("Plots:  output/mainstudy/plots/\n")
cat("RData:  output/mainstudy/N2N_results.RData\n")