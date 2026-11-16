suppressWarnings(library(visNetwork))
library(htmlwidgets)
library(webshot2)
library(tidyverse)
library(PECAN2)

# Output folder for network images
dir.create("output/mainstudy/networks", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Scale validation helper
# Checks that PECAN node values are 1–10 and edge widths are 1–5.
# Prints a warning with the offending rows if anything is out of range.
# Returns TRUE if everything is clean, FALSE if any violations were found.
# =============================================================================
check_pecan_scales <- function(vp, nodes, edges) {
  clean <- TRUE

  bad_values <- nodes[!is.na(nodes$value) & (nodes$value < 1 | nodes$value > 10), ]
  if (nrow(bad_values) > 0) {
    cat("  [SCALE WARNING]", vp, "— node values outside 1–10:\n")
    print(bad_values[, c("id", "label", "value")])
    clean <- FALSE
  }

  bad_widths <- edges[!is.na(edges$width) & (edges$width < 1 | edges$width > 5), ]
  if (nrow(bad_widths) > 0) {
    cat("  [SCALE WARNING]", vp, "— edge widths outside 1–5:\n")
    print(bad_widths[, c("from", "to", "width")])
    clean <- FALSE
  }

  if (clean) cat("  [SCALE OK]", vp, "— all node values and edge widths in range.\n")
  clean
}

# =============================================================================
# Detect which VPs have nar files and which have PECAN files
# =============================================================================
nar_vps <- list.files(
  "data/rawdata/mainstudy/nar",
  pattern = "_nar_dt_nodes\\.csv$"
) |> str_extract("VP\\d+") |> sort()

# Exclude the incomplete VP35 PECAN edge file; only include standard edge files
pecan_edge_vps <- list.files(
  "data/rawdata/mainstudy/PECAN",
  pattern = "^VP\\d+_PECAN_edges\\.csv$"        # excludes "*_unvollständig.csv"
) |> str_extract("VP\\d+")

pecan_node_vps <- list.files(
  "data/rawdata/mainstudy/PECAN",
  pattern = "_PECAN_dt_nodes\\.csv$"
) |> str_extract("VP\\d+")

pecan_vps <- intersect(pecan_node_vps, pecan_edge_vps) |> sort()

cat("VPs with nar files  :", paste(nar_vps,   collapse = ", "), "\n")
cat("VPs with PECAN files:", paste(pecan_vps, collapse = ", "), "\n")

failed <- character(0)

# =============================================================================
# PART 1: Narrative networks (nar_ex and nar_im)
# Loops over all VPs that have nar data files.
#   nar_ex — only nodes/edges where ex.im == "e" (explicitly mentioned)
#   nar_im — all nodes and edges (explicit + implicit)
# =============================================================================

cat("\n--- Narrative networks ---\n")

for (vp in nar_vps) {

  cat("\n[", vp, "] Starting nar networks...\n", sep = "")

  tryCatch({

    # --- Load raw nar data ---
    nar_nodes_all <- read.csv(
      paste0("data/rawdata/mainstudy/nar/", vp, "_nar_dt_nodes.csv"),
      stringsAsFactors = FALSE
    )
    nar_edges_all <- read.csv(
      paste0("data/rawdata/mainstudy/nar/", vp, "_nar_edges.csv"),
      stringsAsFactors = FALSE
    )

    # --- Build nar_ex: explicit nodes and edges only ---
    nar_nodes_ex <- nar_nodes_all |>
      filter(trimws(ex.im) == "e")

    ex_ids <- as.character(nar_nodes_ex$id)

    nar_edges_ex <- nar_edges_all |>
      filter(
        trimws(ex.im) == "e",
        as.character(from) %in% ex_ids,
        as.character(to)   %in% ex_ids
      )

    # --- Build nar_im: all nodes and edges (explicit + implicit) ---
    nar_nodes_im <- nar_nodes_all
    nar_edges_im <- nar_edges_all

    # ---- Visualise nar_ex ------------------------------------------------
    if (nrow(nar_nodes_ex) == 0) {
      cat("[", vp, "] SKIP nar_ex: no explicit nodes found.\n", sep = "")
    } else {
      nar_net_ex <- pecanVis(
        nar_nodes_ex,
        nar_edges_ex,
        nodes_color  = "red",
        edges_color  = "black",
        nodes_label  = "label",
        edges_arrows = "to",
        edges_smooth = TRUE,
        use_seed     = "random",
        c_height     = 750,
        edit         = TRUE,
        export       = TRUE
      ) |>
        visNodes(font = list(color = "black", strokeWidth = 4, strokeColor = "white"))

      tmp_ex <- tempfile(fileext = ".html")
      saveWidget(nar_net_ex, tmp_ex, selfcontained = TRUE)
      webshot(
        tmp_ex,
        file    = paste0("output/mainstudy/networks/", vp, "_nar_ex.png"),
        vwidth  = 1200,
        vheight = 750,
        delay   = 2
      )
      cat("[", vp, "] nar_ex saved.\n", sep = "")
    }

    # ---- Visualise nar_im ------------------------------------------------
    nar_net_im <- pecanVis(
      nar_nodes_im,
      nar_edges_im,
      nodes_color  = "red",
      edges_color  = "black",
      nodes_label  = "label",
      edges_arrows = "to",
      edges_smooth = TRUE,
      use_seed     = "random",
      c_height     = 750,
      edit         = TRUE,
      export       = TRUE
    ) |>
      visNodes(font = list(color = "black", strokeWidth = 4, strokeColor = "white"))

    tmp_im <- tempfile(fileext = ".html")
    saveWidget(nar_net_im, tmp_im, selfcontained = TRUE)
    webshot(
      tmp_im,
      file    = paste0("output/mainstudy/networks/", vp, "_nar_im.png"),
      vwidth  = 1200,
      vheight = 750,
      delay   = 2
    )
    cat("[", vp, "] nar_im saved.\n", sep = "")
    cat("[", vp, "] Done.\n", sep = "")

  }, error = function(e) {
    cat("[", vp, "] ERROR (nar): ", conditionMessage(e), "\n", sep = "")
    failed <<- c(failed, paste0(vp, "_nar"))
  })
}

# =============================================================================
# PART 2: PECAN networks
# Loops over all VPs that have complete PECAN data.
# Edge weights are regulated (reg_by = "width", reg_type = "nodes", reg_value = 1.0),
# matching the convention used in the pilot study visualisation.
# =============================================================================

cat("\n--- PECAN networks ---\n")

for (vp in pecan_vps) {

  cat("\n[", vp, "] Starting PECAN network...\n", sep = "")

  tryCatch({

    pecan_nodes <- read.csv(
      paste0("data/rawdata/mainstudy/PECAN/", vp, "_PECAN_dt_nodes.csv"),
      stringsAsFactors = FALSE
    )
    pecan_edges <- read.csv(
      paste0("data/rawdata/mainstudy/PECAN/", vp, "_PECAN_edges.csv"),
      stringsAsFactors = FALSE
    )

    # --- Scale validation ---
    check_pecan_scales(vp, pecan_nodes, pecan_edges)

    pecan_net <- pecanVis(
      pecan_nodes,
      pecan_edges,
      nodes_size   = "value",
      nodes_color  = "red",
      edges_color  = "black",
      edges_width  = "width",
      nodes_label  = "label",
      edges_arrows = "to",
      edges_smooth = TRUE,
      use_seed     = "random",
      c_width      = 1200,
      c_height     = 750,
      edit         = TRUE,
      export       = TRUE,
      edges_regulation = list(
        reg_by    = "width",
        reg_type  = "nodes",
        reg_value = 1.0
      )
    ) |>
      visNodes(font = list(color = "black", strokeWidth = 4, strokeColor = "white"))

    tmp_pecan <- tempfile(fileext = ".html")
    saveWidget(pecan_net, tmp_pecan, selfcontained = TRUE)
    webshot(
      tmp_pecan,
      file    = paste0("output/mainstudy/networks/", vp, "_pecan.png"),
      vwidth  = 1200,
      vheight = 750,
      delay   = 2
    )
    cat("[", vp, "] PECAN network saved.\n", sep = "")
    cat("[", vp, "] Done.\n", sep = "")

  }, error = function(e) {
    cat("[", vp, "] ERROR (PECAN): ", conditionMessage(e), "\n", sep = "")
    failed <<- c(failed, paste0(vp, "_PECAN"))
  })
}

# =============================================================================
# Summary
# =============================================================================

cat("\n--- Loop complete ---\n")
cat("nar  VPs processed:", length(nar_vps),   "\n")
cat("PECAN VPs processed:", length(pecan_vps), "\n")

if (length(failed) > 0) {
  cat("Failed:", paste(failed, collapse = ", "), "\n")
} else {
  cat("All participants processed successfully.\n")
}

# =============================================================================
# EXAMPLE NETWORKS FOR THESIS: VP24
# =============================================================================
# Three networks saved to output/mainstudy/networks/example/:
#   VP24_nar_ex.png   — explicit nodes only
#   VP24_nar_full.png — full nar network (explicit + implicit)
#   VP24_pecan_reg.png — PECAN with edge regulation (width threshold = 1.0)
#
# Uses VP24_nar_nodes.csv (English labels) rather than VP24_nar_dt_nodes.csv.
# Uses VP24_PECAN_nodes.csv (English labels) rather than VP24_PECAN_dt_nodes.csv.
# Node color: #2E86AB (matches analysis plot theme).
# Label font size: 36 (enlarged for thesis readability).
# =============================================================================

vp         <- "VP24"
node_blue  <- "#2E86AB"   # same blue as N2N_analysis_mainstudy.R plots
label_size <- 36

dir.create("output/mainstudy/networks/example", recursive = TRUE, showWarnings = FALSE)

cat("\n[VP24] Creating thesis example networks...\n")

# --- Load nar data (English labels) ---
nar_nodes_vp24 <- read.csv(
  paste0("data/rawdata/mainstudy/nar/", vp, "_nar_nodes.csv"),
  stringsAsFactors = FALSE
)
nar_edges_vp24 <- read.csv(
  paste0("data/rawdata/mainstudy/nar/", vp, "_nar_edges.csv"),
  stringsAsFactors = FALSE
)

# nar_ex: explicit nodes and edges only
nar_nodes_ex_vp24 <- nar_nodes_vp24 |> filter(trimws(ex.im) == "e")
ex_ids_vp24       <- as.character(nar_nodes_ex_vp24$id)
nar_edges_ex_vp24 <- nar_edges_vp24 |>
  filter(
    trimws(ex.im) == "e",
    as.character(from) %in% ex_ids_vp24,
    as.character(to)   %in% ex_ids_vp24
  )

cat("  nar nodes (full):", nrow(nar_nodes_vp24),    "| edges:", nrow(nar_edges_vp24), "\n")
cat("  nar_ex nodes:    ", nrow(nar_nodes_ex_vp24), "| edges:", nrow(nar_edges_ex_vp24), "\n")

# --- Load PECAN data (English labels from _PECAN_nodes.csv) ---
pecan_nodes_vp24 <- read.csv(
  paste0("data/rawdata/mainstudy/PECAN/", vp, "_PECAN_nodes.csv"),
  stringsAsFactors = FALSE
) |> mutate(label = trimws(label))

pecan_edges_vp24 <- read.csv(
  paste0("data/rawdata/mainstudy/PECAN/", vp, "_PECAN_edges.csv"),
  stringsAsFactors = FALSE
)

cat("  PECAN nodes:     ", nrow(pecan_nodes_vp24), "| edges:", nrow(pecan_edges_vp24), "\n")
check_pecan_scales(vp, pecan_nodes_vp24, pecan_edges_vp24)

# Helper: render a visNetwork widget to PNG
save_net_png <- function(net, filename, vwidth = 1400, vheight = 900) {
  tmp <- tempfile(fileext = ".html")
  saveWidget(net, tmp, selfcontained = TRUE)
  webshot(tmp, file = filename, vwidth = vwidth, vheight = vheight, delay = 2)
  cat("  Saved:", filename, "\n")
}

# Physics helper applied to all three networks:
# - avoidOverlap = 1  : nodes pushed apart until label bounding boxes no longer
#                       collide — the primary fix for label overlap
# - gravitationalConstant: more negative = stronger repulsion between nodes
# - springLength: longer resting length = more breathing room between nodes
# - stabilization iterations: more time for the layout to fully settle
spread_physics <- function(net) {
  visPhysics(
    net,
    solver = "barnesHut",
    barnesHut = list(
      gravitationalConstant = -12000,
      centralGravity        = 0.5,  # higher value keeps isolated nodes near centre
      springLength          = 200,
      avoidOverlap          = 1     # key: prevents label bounding-box overlap
    ),
    stabilization = list(enabled = TRUE, iterations = 400)
  )
}

# Shared font style: black text with white stroke (readable on blue nodes)
node_font <- list(size = label_size, color = "black",
                  strokeWidth = 5, strokeColor = "white")

# --- 1. nar_ex network (explicit nodes only) ---
nar_net_ex_vp24 <- pecanVis(
  nar_nodes_ex_vp24,
  nar_edges_ex_vp24,
  nodes_color  = node_blue,
  edges_color  = "black",
  nodes_label  = "label",
  edges_arrows = "to",
  edges_smooth = TRUE,
  use_seed     = "random",
  c_width      = 1400,
  c_height     = 900,
  edit         = TRUE,
  export       = TRUE
) |>
  visNodes(font = node_font) |>
  spread_physics()

save_net_png(
  nar_net_ex_vp24,
  paste0("output/mainstudy/networks/example/", vp, "_nar_ex.png")
)

# --- 2. nar full network (explicit + implicit) ---
nar_net_full_vp24 <- pecanVis(
  nar_nodes_vp24,
  nar_edges_vp24,
  nodes_color  = node_blue,
  edges_color  = "black",
  nodes_label  = "label",
  edges_arrows = "to",
  edges_smooth = TRUE,
  use_seed     = "random",
  c_width      = 1400,
  c_height     = 900,
  edit         = TRUE,
  export       = TRUE
) |>
  visNodes(font = node_font) |>
  spread_physics()

save_net_png(
  nar_net_full_vp24,
  paste0("output/mainstudy/networks/example/", vp, "_nar_full.png")
)

# --- 3. PECAN regulated network ---
pecan_net_vp24 <- pecanVis(
  pecan_nodes_vp24,
  pecan_edges_vp24,
  nodes_size   = "value",
  nodes_color  = node_blue,
  edges_color  = "black",
  edges_width  = "width",
  nodes_label  = "label",
  edges_arrows = "to",
  edges_smooth = TRUE,
  use_seed     = "random",
  c_width      = 1400,
  c_height     = 900,
  edit         = TRUE,
  export       = TRUE,
  edges_regulation = list(
    reg_by    = "width",
    reg_type  = "nodes",
    reg_value = 1.0
  )
) |>
  visNodes(font = node_font) |>
  spread_physics()

save_net_png(
  pecan_net_vp24,
  paste0("output/mainstudy/networks/example/", vp, "_pecan_reg.png")
)

cat("[VP24] All three example networks complete.\n")