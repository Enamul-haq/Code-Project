# ============================================================
# Rheumatoid Arthritis: Immune Infiltration (xCell) + Hub Gene Correlation
# Datasets (CSV ready):
#   Discovery : GSE93272  (whole blood)
#   Validation: GSE45291  (whole blood; RA only; SLE বাদ)
#
# OUTPUTS (TIFF 600 DPI):
#   Fig1A : Heatmap – Hub genes vs Immune cells (Discovery)
#   Fig1B : Heatmap – Hub genes vs Immune cells (Validation)
#   Fig2  : Bubble plot – Both cohorts combined
#   Fig4  : Violin panels – Top 18 immune cells (each dataset)
#   FigBox: Reference-style boxplot – 18 immune cells (each dataset)
#
# TABLES (CSV):
#   Table1_Immune_Wilcoxon.csv
#   Table2_HubGene_Immune_Correlation.csv
# ============================================================

rm(list = ls())
set.seed(123)

# ─────────────────────────────────────────────────────────
# 0)  USER SETTINGS  ← edit only this block
# ─────────────────────────────────────────────────────────
setwd("E:/2 Rheumatoid Arthritis/7 Immune infiltration analysis")  # <<< CHANGE to your folder

files <- c(
  "cleaned_modified_final_data_GSE93272.csv",  # Discovery
  "cleaned_modified_final_data_GSE45291.csv"   # Validation
)

ds_names <- c("GSE93272", "GSE45291")
ds_labels <- c(
  "Immune infiltration (xCell): GSE93272 (RA)",
  "Immune infiltration (xCell): GSE45291 (RA)"
)

# RA Key genes (KGs): Down + Up (your list)
HUB_GENES <- c(
  "TNFAIP3","JUN","MYC","CTNNB1","RHOA",   # down
  "IL1B","TLR4","CD86","FCGR3B","CD8A"     # up
)

OUTDIR <- "RA_ImmuneKG_xCell_Publication"

# Group label (output plots/tables will show these)
CTRL_LABEL <- "Control"
CASE_LABEL <- "RA"

# Sample column naming rule (edit only if your columns differ)
# Default expects: control..., disease...
CTRL_PATTERN <- "^control"
CASE_PATTERN <- "^disease"

# Colors (used in some plots)
COL_CTRL <- "#E07070"  # salmon-red
COL_CASE <- "#4DBDBD"  # teal

# Heatmap colors for gene strip (auto recycle if needed)
GENE_COLOURS <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF","#00BA38",
                  "#619CFF","#F564E3","#A3A500","#00A9FF","#B79F00",
                  "#00A087","#4DBBD5","#E64B35","#3C5488")

# ─────────────────────────────────────────────────────────
# 1)  PACKAGES
# ─────────────────────────────────────────────────────────
# =========================
# Smart install (only missing)
# =========================

install_if_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing (CRAN): ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  } else {
    message("All CRAN packages already installed.")
  }
}

install_if_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing (Bioconductor): ", paste(missing, collapse = ", "))
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    message("All Bioconductor packages already installed.")
  }
}

install_xcell_if_missing <- function() {
  if (!requireNamespace("xCell", quietly = TRUE)) {
    # needs remotes
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    message("Installing xCell from GitHub (dviraran/xCell) ...")
    remotes::install_github("dviraran/xCell", upgrade = "never")
  } else {
    message("xCell already installed.")
  }
}

# ---- 1) CRAN packages ----
cran_pkgs <- c(
  "data.table","dplyr","tibble","stringr","tidyr",
  "ggplot2","ggpubr","ragg","patchwork","forcats","scales",
  "circlize"
)
install_if_missing_cran(cran_pkgs)

# ---- 2) Bioconductor packages ----
bioc_pkgs <- c("ComplexHeatmap")
install_if_missing_bioc(bioc_pkgs)

# ---- 3) GitHub package ----
install_xcell_if_missing()

# =========================
# Now load packages safely
# =========================
suppressPackageStartupMessages({
  library(data.table);   library(dplyr);     library(tibble)
  library(stringr);      library(tidyr);     library(ggplot2)
  library(ggpubr);       library(ragg);      library(patchwork)
  library(xCell);        library(circlize);  library(ComplexHeatmap)
  library(grid);         library(forcats);   library(scales)
})

# ─────────────────────────────────────────────────────────
# 2)  OUTPUT FOLDERS
# ─────────────────────────────────────────────────────────
fig_dir <- file.path(OUTDIR, "figures_tiff")
tab_dir <- file.path(OUTDIR, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────
# 3)  HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(face = "bold", size = base_size - 1),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      strip.background = element_rect(fill = "grey88", colour = "grey50"),
      plot.title       = element_text(face = "bold", hjust = 0.5,
                                      size = base_size + 2),
      legend.title     = element_text(face = "bold")
    )
}

save_tiff <- function(plt, filename, w = 14, h = 8, dpi = 600) {
  ggsave(file.path(fig_dir, filename), plot = plt,
         device = "tiff", dpi = dpi,
         width = w, height = h, units = "in",
         compression = "lzw")
  message("  Saved: ", filename)
}

detect_gene_col <- function(df) {
  cands <- c("gene_symbol","GeneSymbol","Gene","SYMBOL","symbol","gene")
  hit   <- intersect(cands, colnames(df))
  if (!length(hit)) stop("No gene-symbol column found. Need one of: gene_symbol/Gene/SYMBOL etc.")
  hit[1]
}

make_expr_and_meta <- function(file,
                               ctrl_label = CTRL_LABEL,
                               case_label = CASE_LABEL,
                               ctrl_pat = CTRL_PATTERN,
                               case_pat = CASE_PATTERN) {
  
  df   <- fread(file, data.table = FALSE)
  gcol <- detect_gene_col(df)
  colnames(df)[colnames(df) == gcol] <- "gene_symbol"
  
  df <- df %>% filter(!is.na(gene_symbol), gene_symbol != "")
  
  sc <- setdiff(colnames(df), "gene_symbol")
  
  if (!any(str_detect(sc, ctrl_pat)) || !any(str_detect(sc, case_pat))) {
    stop(paste0(
      "Sample columns must match patterns.\n",
      "CTRL_PATTERN=", ctrl_pat, "  CASE_PATTERN=", case_pat, "\n",
      "Example required: control_1... and disease_1...\n",
      "Your sample columns are: ", paste(head(sc, 10), collapse = ", "), " ..."
    ))
  }
  
  df[sc] <- lapply(df[sc], function(x) as.numeric(as.character(x)))
  
  # collapse duplicate genes
  df <- df %>%
    group_by(gene_symbol) %>%
    summarise(across(all_of(sc), ~ median(.x, na.rm = TRUE)), .groups = "drop")
  
  expr <- df %>% column_to_rownames("gene_symbol") %>% as.matrix()
  mode(expr) <- "numeric"
  
  grp <- ifelse(str_detect(colnames(expr), ctrl_pat), ctrl_label,
                ifelse(str_detect(colnames(expr), case_pat), case_label, NA))
  
  meta <- data.frame(
    sample = colnames(expr),
    group  = factor(grp, levels = c(ctrl_label, case_label)),
    stringsAsFactors = FALSE
  )
  
  list(expr = expr, meta = meta)
}

wilcox_cells <- function(sm, meta, ctrl_label = CTRL_LABEL, case_label = CASE_LABEL) {
  lapply(rownames(sm), function(ct) {
    x  <- as.numeric(sm[ct, meta$group == ctrl_label])
    y  <- as.numeric(sm[ct, meta$group == case_label])
    wt <- wilcox.test(x, y, exact = FALSE)
    
    data.frame(
      cell        = ct,
      p           = wt$p.value,
      med_control = median(x, na.rm = TRUE),
      med_case    = median(y, na.rm = TRUE),
      delta_median= median(y, na.rm = TRUE) - median(x, na.rm = TRUE)
    )
  }) %>% bind_rows() %>%
    mutate(FDR = p.adjust(p, "BH")) %>%
    arrange(FDR, p)
}

corr_kg_immune <- function(expr, immune, genes, cells) {
  gu <- intersect(genes, rownames(expr))
  cu <- intersect(cells, rownames(immune))
  if (!length(gu)) stop("Hub genes not found in expression matrix. Check gene symbols.")
  if (!length(cu)) stop("Immune cells not found in xCell output.")
  
  expand.grid(gene = gu, cell = cu, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      rho = suppressWarnings(cor(as.numeric(expr[gene, ]),
                                 as.numeric(immune[cell, ]),
                                 method = "spearman")),
      p   = suppressWarnings(cor.test(as.numeric(expr[gene, ]),
                                      as.numeric(immune[cell, ]),
                                      method = "spearman")$p.value)
    ) %>%
    ungroup() %>%
    mutate(FDR = p.adjust(p, "BH"))
}

make_heatmap_like_ref <- function(cor_df, dataset_label, filename,
                                  gene_order = NULL, col_order = NULL,
                                  panel_tag = "A",
                                  lim = 0.6,
                                  row_k = 3, col_k = 3) {
  
  mat <- cor_df %>%
    dplyr::select(gene, cell, rho) %>%
    tidyr::pivot_wider(names_from = cell, values_from = rho) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  if (!is.null(col_order)) {
    keepC <- intersect(col_order, colnames(mat))
    mat <- mat[, keepC, drop = FALSE]
  }
  if (!is.null(gene_order)) {
    keepR <- intersect(gene_order, rownames(mat))
    mat <- mat[keepR, , drop = FALSE]
  }
  
  gl <- rownames(mat)
  gcols <- setNames(GENE_COLOURS[seq_len(length(gl))], gl)
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Gene = gl,
    col = list(Gene = gcols),
    show_annotation_name = FALSE,
    annotation_width = grid::unit(4, "mm")
  )
  
  col_fun <- circlize::colorRamp2(
    c(-lim, 0, lim),
    c("#0000FF", "white", "#FF0000")
  )
  
  ragg::agg_tiff(
    filename = file.path(fig_dir, filename),
    width = 12, height = 7.5, units = "in",
    res = 600, compression = "lzw"
  )
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = NULL,
    col  = col_fun,
    
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    
    row_split = row_k,
    column_split = col_k,
    row_gap = grid::unit(2.5, "mm"),
    column_gap = grid::unit(2.5, "mm"),
    
    left_annotation = row_ha,
    
    row_names_gp = grid::gpar(fontface = "plain", fontsize = 13),
    column_names_gp = grid::gpar(fontface = "plain", fontsize = 12),
    column_names_rot = 90,
    column_names_side = "bottom",
    
    rect_gp = grid::gpar(col = "grey70", lwd = 1.3),
    border = TRUE,
    
    row_dend_width = grid::unit(18, "mm"),
    column_dend_height = grid::unit(14, "mm"),
    
    heatmap_legend_param = list(
      title = "rho",
      title_position = "topcenter",
      title_gp  = grid::gpar(fontface = "plain", fontsize = 12),
      at = c(-lim, -0.4, -0.2, 0, 0.2, 0.4, lim),
      labels_gp = grid::gpar(fontsize = 12),
      legend_height = grid::unit(55, "mm"),
      legend_width  = grid::unit(10, "mm"),
      tick_length = grid::unit(0, "mm")
    ),
    
    column_title = paste0("(", panel_tag, ") Correlation: Hub Genes vs Immune Cells (", dataset_label, ")"),
    column_title_gp = grid::gpar(fontface = "bold", fontsize = 16)
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE
  )
  
  dev.off()
  message("  Saved: ", filename)
}

plot_all_cells <- function(immune, meta, dataset_label, cells_use,
                           ctrl_label = CTRL_LABEL, case_label = CASE_LABEL) {
  
  df <- as.data.frame(t(immune[cells_use, , drop = FALSE]))
  df$sample <- rownames(df)
  
  df_long <- df %>%
    left_join(meta, by = "sample") %>%
    pivot_longer(cols = all_of(cells_use),
                 names_to = "Cell",
                 values_to = "Score") %>%
    mutate(
      Cell  = fct_relevel(Cell, cells_use),
      group = fct_relevel(group, c(ctrl_label, case_label))
    )
  
  ggplot(df_long, aes(x = group, y = Score, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.85) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95) +
    stat_compare_means(method = "wilcox.test",
                       label = "p.format",
                       size = 3) +
    facet_wrap(~Cell, scales = "free_y", ncol = 6) +
    labs(title = dataset_label, x = NULL, y = "xCell score") +
    theme_pub(11) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

plot_immune_box_refstyle <- function(immune, meta, dataset_label, cells_use,
                                     ctrl_label = CTRL_LABEL, case_label = CASE_LABEL,
                                     ylab_text = "xCell score") {
  
  COL_CTRL_REF <- "#0072B2"  # blue
  COL_DIS_REF  <- "#D55E00"  # orange/gold
  
  df <- as.data.frame(t(immune[cells_use, , drop = FALSE]))
  df$sample <- rownames(df)
  
  df_long <- df %>%
    left_join(meta, by = "sample") %>%
    pivot_longer(cols = all_of(cells_use), names_to = "Cell", values_to = "Score") %>%
    mutate(
      Cell  = fct_relevel(Cell, cells_use),
      group = fct_relevel(group, c(ctrl_label, case_label))
    )
  
  ggplot(df_long, aes(x = Cell, y = Score, fill = group, colour = group)) +
    geom_boxplot(
      width = 0.55,
      position = position_dodge(width = 0.75),
      outlier.shape = NA,
      alpha = 1
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.75),
      size = 1.8, alpha = 0.95, stroke = 0
    ) +
    scale_fill_manual(values = setNames(c(COL_CTRL_REF, COL_DIS_REF), c(ctrl_label, case_label)), name = "Group") +
    scale_colour_manual(values = setNames(c(COL_CTRL_REF, COL_DIS_REF), c(ctrl_label, case_label)), guide = "none") +
    labs(title = dataset_label, x = NULL, y = ylab_text) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.title.y = element_text(face = "plain"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(6, "mm")
    )
}

# ─────────────────────────────────────────────────────────
# 4)  LOAD DATA  +  RUN xCell
# ─────────────────────────────────────────────────────────
message("Loading data …")
d1 <- make_expr_and_meta(files[1])
d2 <- make_expr_and_meta(files[2])

message("Running xCell: ", ds_names[1])
imm1 <- xCellAnalysis(d1$expr)
imm1 <- imm1[, d1$meta$sample, drop = FALSE]

message("Running xCell: ", ds_names[2])
imm2 <- xCellAnalysis(d2$expr)
imm2 <- imm2[, d2$meta$sample, drop = FALSE]

# ─────────────────────────────────────────────────────────
# 5)  IMMUNE PANEL SELECTION
# ─────────────────────────────────────────────────────────
stats1 <- wilcox_cells(imm1, d1$meta)
stats2 <- wilcox_cells(imm2, d2$meta)

panel_cells <- stats1 %>%
  filter(cell %in% rownames(imm2)) %>%
  pull(cell) %>%
  head(18)

top10_cells <- head(panel_cells, 10)

message("Panel cells (", length(panel_cells), "): ",
        paste(panel_cells, collapse = ", "))

# ─────────────────────────────────────────────────────────
# 6)  TABLE 1 – Wilcoxon
# ─────────────────────────────────────────────────────────
tab1 <- bind_rows(
  stats1 %>% mutate(dataset = ds_names[1]),
  stats2 %>% mutate(dataset = ds_names[2])
) %>%
  filter(cell %in% panel_cells) %>%
  select(dataset, cell, med_control, med_case, delta_median, p, FDR) %>%
  arrange(cell, dataset)

write.csv(tab1, file.path(tab_dir, "Table1_Immune_Wilcoxon.csv"), row.names = FALSE)
message("  Saved: Table1_Immune_Wilcoxon.csv")

# ─────────────────────────────────────────────────────────
# 7)  TABLE 2 – Hub gene × Immune correlation
# ─────────────────────────────────────────────────────────
cor1 <- corr_kg_immune(d1$expr, imm1, HUB_GENES, panel_cells) %>% mutate(dataset = ds_names[1])
cor2 <- corr_kg_immune(d2$expr, imm2, HUB_GENES, panel_cells) %>% mutate(dataset = ds_names[2])

tab2 <- bind_rows(cor1, cor2) %>% arrange(dataset, gene, cell)
write.csv(tab2, file.path(tab_dir, "Table2_HubGene_Immune_Correlation.csv"), row.names = FALSE)
message("  Saved: Table2_HubGene_Immune_Correlation.csv")

# ─────────────────────────────────────────────────────────
# 8)  FIG 1A & 1B – ComplexHeatmap (publication style)
# ─────────────────────────────────────────────────────────
make_heatmap_like_ref(
  cor1, ds_names[1], "Fig1A_Heatmap_GSE93272.tiff",
  gene_order = HUB_GENES, col_order = panel_cells,
  panel_tag = "A", row_k = 3, col_k = 3
)

make_heatmap_like_ref(
  cor2, ds_names[2], "Fig1B_Heatmap_GSE45291.tiff",
  gene_order = HUB_GENES, col_order = panel_cells,
  panel_tag = "B", row_k = 3, col_k = 3
)

# ─────────────────────────────────────────────────────────
# 9)  FIG 2 – Bubble plot (both cohorts)
# ─────────────────────────────────────────────────────────
lim <- 0.6
col_breaks <- seq(-lim, lim, by = 0.2)
fmt_dec <- function(x) formatC(x, format = "f", digits = 1)

bubble_df <- tab2 %>%
  mutate(
    abs_rho = abs(rho),
    sig     = ifelse(FDR < 0.05, "FDR<0.05", "ns"),
    dataset = factor(dataset, levels = ds_names),
    cell    = factor(cell, levels = panel_cells),
    gene    = factor(gene, levels = rev(HUB_GENES))
  )

p_bubble <- ggplot(
  bubble_df,
  aes(x = cell, y = gene, size = abs_rho, colour = rho, shape = dataset, alpha = sig)
) +
  geom_point(stroke = 0.6) +
  scale_alpha_manual(values = c("FDR<0.05" = 1.0, "ns" = 0.12), guide = "none") +
  scale_colour_gradient2(
    low  = "#0000FF", mid = "white", high = "#FF0000",
    midpoint = 0, limits = c(-lim, lim),
    breaks = col_breaks, labels = fmt_dec,
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(
    range = c(2.2, 8.2),
    breaks = c(0.2, 0.4, 0.6),
    labels = fmt_dec,
    name = "|rho|"
  ) +
  scale_shape_manual(values = c(16, 17), name = "Cohort") +
  labs(
    title = "Hub Gene–Immune Cell Correlation (xCell; Spearman)",
    x = "Immune cell type",
    y = "Key gene"
  ) +
  theme_pub(12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    legend.box = "vertical",
    legend.spacing.y = unit(4, "mm"),
    legend.key.height = unit(7, "mm"),
    legend.key.width  = unit(7, "mm")
  ) +
  guides(
    colour = guide_colourbar(
      order = 1,
      barheight = unit(45, "mm"),
      barwidth  = unit(8, "mm"),
      ticks = FALSE,
      frame.colour = "white"
    ),
    shape = guide_legend(order = 2, override.aes = list(alpha = 1, size = 3)),
    size  = guide_legend(order = 3, override.aes = list(alpha = 1, shape = 16))
  )

save_tiff(p_bubble, "Fig2_Bubble_BothCohorts.tiff", w = 15, h = 7)

# ─────────────────────────────────────────────────────────
# 10)  FIG 4 – Violin panels (18 cells) + reference boxplots
# ─────────────────────────────────────────────────────────
p_all_A <- plot_all_cells(imm1, d1$meta, paste0(ds_labels[1]), panel_cells)
save_tiff(p_all_A, paste0("Fig4_", ds_names[1], "_All18Cells.tiff"), w = 16, h = 9)

p_all_B <- plot_all_cells(imm2, d2$meta, paste0(ds_labels[2]), panel_cells)
save_tiff(p_all_B, paste0("Fig4_", ds_names[2], "_All18Cells.tiff"), w = 16, h = 9)

p_box_A <- plot_immune_box_refstyle(imm1, d1$meta, paste0(ds_names[1], " (xCell)"), panel_cells)
save_tiff(p_box_A, paste0("FigBox_All18Cells_", ds_names[1], ".tiff"), w = 16, h = 7)

p_box_B <- plot_immune_box_refstyle(imm2, d2$meta, paste0(ds_names[2], " (xCell)"), panel_cells)
save_tiff(p_box_B, paste0("FigBox_All18Cells_", ds_names[2], ".tiff"), w = 16, h = 7)

message("\nDONE ✅ Outputs saved under: ", OUTDIR)
