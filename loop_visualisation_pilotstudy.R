suppressWarnings(library(visNetwork))
library(htmlwidgets)
library(webshot2)
library(PECAN2)

# Output folder for network images
dir.create("output/pilotstudy/networks", recursive = TRUE, showWarnings = FALSE)

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

# Replace ae/oe/ue substitutions with proper German umlauts
fix_umlauts <- function(x) {
  x |>
    gsub("UE", "Ü", x = _, fixed = TRUE) |>
    gsub("Ue", "Ü", x = _, fixed = TRUE) |>
    gsub("OE", "Ö", x = _, fixed = TRUE) |>
    gsub("Oe", "Ö", x = _, fixed = TRUE) |>
    gsub("AE", "Ä", x = _, fixed = TRUE) |>
    gsub("Ae", "Ä", x = _, fixed = TRUE) |>
    gsub("ue", "ü", x = _, fixed = TRUE) |>
    gsub("oe", "ö", x = _, fixed = TRUE) |>
    gsub("ae", "ä", x = _, fixed = TRUE)
}

failed <- character(0)

for (vp in sprintf("VP%02d", 1:22)) {
  
  cat("\n[", vp, "] Starting...\n", sep = "")
  
  tryCatch({
    
    # --- Load data ---
    pecan_nodes <- read.csv(paste0("data/rawdata/pilotstudy/PECAN/", vp, "_PECAN_dt_nodes.csv"))
    pecan_edges <- read.csv(paste0("data/rawdata/pilotstudy/PECAN/", vp, "_PECAN_edges.csv"))
    nar_nodes   <- read.csv(paste0("data/rawdata/pilotstudy/nar/",   vp, "_nar_dt_nodes.csv"))
    nar_edges   <- read.csv(paste0("data/rawdata/pilotstudy/nar/",   vp, "_nar_edges.csv"))
    cat("[", vp, "] Data loaded.\n", sep = "")
    
    # --- Fix umlauts in labels ---
    pecan_nodes$label <- fix_umlauts(pecan_nodes$label)
    nar_nodes$label   <- fix_umlauts(nar_nodes$label)

    # --- Scale validation ---
    check_pecan_scales(vp, pecan_nodes, pecan_edges)

    # --- Narrative network ---
    nar_net <- pecanVis(nar_nodes,
                        nar_edges,
                        nodes_color  = "red",
                        edges_color  = "black",
                        nodes_label  = "label",
                        edges_arrows = "to",
                        edges_smooth = TRUE,
                        use_seed     = "random",
                        c_height     = 750,
                        edit         = TRUE,
                        export       = TRUE) |>
      visNodes(font = list(color = "black", strokeWidth = 4, strokeColor = "white"))
    
    tmp_nar <- tempfile(fileext = ".html")
    saveWidget(nar_net, tmp_nar, selfcontained = TRUE)
    webshot(tmp_nar,
            file    = paste0("output/pilotstudy/networks/", vp, "_nar.png"),
            vwidth  = 1200,
            vheight = 750)
    cat("[", vp, "] Narrative network saved as PNG.\n", sep = "")
    
    # --- PECAN network ---
    pecan_net <- pecanVis(pecan_nodes,
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
                          )) |>
      visNodes(font = list(color = "black", strokeWidth = 4, strokeColor = "white"))
    
    tmp_pecan <- tempfile(fileext = ".html")
    saveWidget(pecan_net, tmp_pecan, selfcontained = TRUE)
    webshot(tmp_pecan,
            file    = paste0("output/pilotstudy/networks/", vp, "_pecan.png"),
            vwidth  = 1200,
            vheight = 750)
    cat("[", vp, "] PECAN network saved as PNG.\n", sep = "")
    
    cat("[", vp, "] Done.\n", sep = "")
    
  }, error = function(e) {
    cat("[", vp, "] ERROR: ", conditionMessage(e), "\n", sep = "")
    failed <<- c(failed, vp)
  })
}

cat("\n--- Loop complete ---\n")
if (length(failed) > 0) {
  cat("Failed participants:", paste(failed, collapse = ", "), "\n")
} else {
  cat("All 22 participants processed successfully.\n")
}