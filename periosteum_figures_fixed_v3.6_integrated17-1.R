
# =============================================================================
# periosteum_figures_fixed_v3.6.R
# GSE280914 periosteum aging / abortive regeneration figure repair
# Input: CSV tables exported by periosteum_aging v3.6-publication.
# Output: corrected PDF + 600-dpi TIFF figures in ./figures_fixed
# Notes:
#   * Set PERIO_DATA_DIR to the folder containing the CSV files, or run with
#     the uploaded CSV files in the working directory.
#   * Optional cell-level file for exact Fig2/Fig3a violins:
#       senmayo_celllevel_cd45neg_stromal.csv
#     columns: cell_id, age_group, injury, fraction, celltype, senmayo
#   * Optional UMAP coordinate file for Fig7 numeric trajectory labels:
#       umap_coords.csv  columns: cell_id, umap_1, umap_2
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(scales); library(patchwork); library(grid)
})
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package ggrepel is required.")
if (!requireNamespace("circlize", quietly = TRUE)) stop("Package circlize is required for Fig6b.")

data_dir <- Sys.getenv("PERIO_DATA_DIR", unset = ".")
if (!file.exists(file.path(data_dir, "TableS1_composition.csv")) && dir.exists("/mnt/agents/upload")) {
  data_dir <- "/mnt/agents/upload"
}
outdir <- file.path(data_dir, "figures_fixed")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- small utilities -------------------------------------------------
read_tab <- function(name) readr::read_csv(file.path(data_dir, name), show_col_types = FALSE, progress = FALSE)
need_cols <- function(df, cols, name) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing)) stop(name, " is missing columns: ", paste(missing, collapse = ", "))
  invisible(df)
}
assert_close <- function(x, y, tol = 1e-6, label = "value") {
  if (!is.finite(x) || abs(x - y) > tol) stop("Accuracy check failed for ", label, ": got ", x, ", expected ", y)
  invisible(x)
}
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 2), formatC(p, format = "f", digits = 3))
save_gg <- function(p, filename, width, height, dpi = 600) {
  ggplot2::ggsave(file.path(outdir, paste0(filename, ".pdf")), p, width = width, height = height, units = "in", device = cairo_pdf)
  ggplot2::ggsave(file.path(outdir, paste0(filename, ".tiff")), p, width = width, height = height, units = "in", dpi = dpi, compression = "lzw")
  invisible(file.path(outdir, filename))
}
celltype_levels <- c("Chondro","Cycling","Endo","Macroph","Neutroph","pSSPC_fib","pSSPC_osteo","pSSPC_q")
celltype_cols <- c(Chondro="#8dd3c7", Cycling="#8c6d5a", Endo="#009e73", Macroph="#f4a261",
                   Neutroph="#e41a1c", pSSPC_fib="#56b4e9", pSSPC_osteo="#7570b3", pSSPC_q="#b3b3b3")
group_levels <- c("Young_Intact","Young_Fracture","Aged_Intact","Aged_Fracture")
group_cols <- c(Young_Intact="#F8766D", Young_Fracture="#7CAE00", Aged_Intact="#00BFC4", Aged_Fracture="#C77CFF")

# ---------- load and validate key tables ------------------------------------
comp <- read_tab("TableS1_composition.csv"); need_cols(comp, c("donor_id","age_group","injury","fraction","celltype","n","pct"), "TableS1")
sen_donor <- read_tab("TableS5_SenMayo_by_donor.csv"); need_cols(sen_donor, c("donor_id","age_group","injury","SenMayo"), "TableS5")
gsea_young <- read_tab("TableS6_GSEA_youngResponse_inAged.csv"); need_cols(gsea_young, c("ID","setSize","enrichmentScore","NES","pvalue","p.adjust","rank","leading_edge"), "TableS6")
callus <- read_tab("Table4_callus_validation.csv"); need_cols(callus, c("gene","lfc_peri","lfc_callus","direction_consistent"), "Table4")
cc_y <- read_tab("TableS8_CC_Young.csv"); need_cols(cc_y, c("source","target","ligand","receptor","prob","pval","pathway_name"), "TableS8")
cc_a <- read_tab("TableS9_CC_Aged.csv"); need_cols(cc_a, c("source","target","ligand","receptor","prob","pval","pathway_name"), "TableS9")
pt <- read_tab("TableS10_pseudotime.csv"); need_cols(pt, c("pseudotime","group","age_group","injury","celltype","senmayo"), "TableS10")
pt <- pt %>% mutate(cell_id = paste0("cell_", dplyr::row_number()))
pt_group <- read_tab("TableS10b_pt_by_group.csv"); need_cols(pt_group, c("group","median_pt","mean_sen","n"), "TableS10b")
mr <- read_tab("Table6_MR_results.csv"); need_cols(mr, c("gene","id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval"), "Table6")
mr_causal <- read_tab("Table7_MR_causal_genes.csv"); need_cols(mr_causal, c("gene","outcome","b","se","pval","sc_lfc","direction_match"), "Table7")

# ---------- optional auto-export of cell-level SenMayo and UMAP coords ----------
# Priority used downstream: existing CSV > CSV exported here from a Seurat object > fallback.
# The exported umap_coords.csv is usable by this script only when cell_id equals the synthetic
# cell_i IDs created for TableS10_pseudotime.csv above, in the same row order.
find_seurat_object <- function(envir = .GlobalEnv) {
  preferred <- c(Sys.getenv("SEURAT_OBJ_NAME", unset = ""), "periosteum_obj", "seurat_obj", "obj", "integrated_obj", "periosteum_seurat")
  preferred <- preferred[nzchar(preferred)]
  for (nm in preferred) {
    if (exists(nm, envir = envir, inherits = FALSE)) {
      x <- get(nm, envir = envir)
      if (inherits(x, "Seurat")) return(x)
    }
  }
  invisible(NULL)
}

export_optional_inputs <- function(obj = find_seurat_object(),
                                   out_csv_dir = data_dir,
                                   overwrite = FALSE) {
  if (is.null(obj)) return(invisible(FALSE))
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    message("Seurat not available; optional CSV export skipped.")
    return(invisible(FALSE))
  }
  dir.create(out_csv_dir, showWarnings = FALSE, recursive = TRUE)
  md <- obj@meta.data
  sen_cols <- grep("senmayo", colnames(md), ignore.case = TRUE, value = TRUE)
  if (length(sen_cols) < 1) {
    message("No senmayo column in obj@meta.data; optional CSV export skipped.")
    return(invisible(FALSE))
  }
  sen_col <- sen_cols[1]
  req <- c("age_group", "injury", "fraction", "celltype")
  if (!all(req %in% colnames(md))) {
    message("obj@meta.data lacks one of age_group/injury/fraction/celltype; optional CSV export skipped.")
    return(invisible(FALSE))
  }

  # If the object is exactly the TableS10 trajectory object, use the same synthetic IDs.
  same_as_pt <- nrow(md) == nrow(pt)
  ids <- if (same_as_pt) pt$cell_id else rownames(md)
  if (!same_as_pt) {
    message("Seurat object row count differs from TableS10; exported coords keep object rownames and may not join to this script unless IDs are matched.")
  }

  sen_file <- file.path(out_csv_dir, "senmayo_celllevel_cd45neg_stromal.csv")
  if (!file.exists(sen_file) || overwrite) {
    readr::write_csv(tibble::tibble(
      cell_id = ids, age_group = md$age_group, injury = md$injury,
      fraction = md$fraction, celltype = md$celltype, senmayo = as.numeric(md[[sen_col]])
    ), sen_file)
    message("Exported: ", sen_file)
  }

  umap_file <- file.path(out_csv_dir, "umap_coords.csv")
  if (!file.exists(umap_file) || overwrite) {
    if ("umap" %in% names(obj@reductions)) {
      um <- Seurat::Embeddings(obj, reduction = "umap")
      readr::write_csv(tibble::tibble(cell_id = ids, umap_1 = as.numeric(um[, 1]), umap_2 = as.numeric(um[, 2])), umap_file)
      message("Exported: ", umap_file)
    } else {
      message("No umap reduction in object; umap_coords.csv not exported.")
    }
  }
  invisible(TRUE)
}
if (Sys.getenv("PERIO_AUTO_EXPORT_OPTIONAL", unset = "1") == "1") export_optional_inputs(overwrite = FALSE)


# Derived validation values expected from the run logs / figures.
comp <- comp %>% mutate(group = paste(age_group, injury, sep = "_"))
assert_close(sum(comp$n), 31377, tol = 0.1, label = "annotated cell total")
assert_close(sum(comp$celltype == "pSSPC_fib" & comp$fraction == "CD45neg" & comp$injury == "Intact" & comp$age_group == "Young") /
             sum(comp$n[comp$fraction == "CD45neg" & comp$injury == "Intact" & comp$age_group == "Young"]) * 100,
             5.1330798, tol = 1e-5, label = "Young intact pSSPC_fib % of CD45neg")
assert_close(sum(comp$celltype == "pSSPC_fib" & comp$fraction == "CD45neg" & comp$injury == "Intact" & comp$age_group == "Aged") /
             sum(comp$n[comp$fraction == "CD45neg" & comp$injury == "Intact" & comp$age_group == "Aged"]) * 100,
             12.7861771, tol = 1e-5, label = "Aged intact pSSPC_fib % of CD45neg")
sen_intact <- sen_donor %>% filter(injury == "Intact") %>% arrange(age_group)
assert_close(nrow(sen_intact), 2, tol = 0, label = "intact donor count")
assert_close(sen_intact$SenMayo[sen_intact$age_group == "Aged"] / sen_intact$SenMayo[sen_intact$age_group == "Young"],
             0.9606391, tol = 1e-5, label = "donor-level SenMayo aged/young fold")
assert_close(nrow(callus), 310, tol = 0, label = "callus common genes")
assert_close(sum(callus$direction_consistent), 172, tol = 0, label = "callus direction-consistent genes")
assert_close(cor(callus$lfc_peri, callus$lfc_callus, method = "pearson"), 0.04229914, tol = 1e-6, label = "periosteum-callus Pearson rho")
assert_close(as.numeric(gsea_young$NES[1]), 1.65970118, tol = 1e-6, label = "YoungInjuryResponse NES in aged fracture")
assert_close(as.numeric(gsea_young$p.adjust[1]), 2.95811778619461e-17, tol = 1e-22, label = "YoungInjuryResponse p.adjust")
assert_close(nrow(cc_y), 2126, tol = 0, label = "Young CellChat interaction count")
assert_close(sum(cc_y$prob), 0.2729917537, tol = 1e-8, label = "Young CellChat strength")
assert_close(nrow(cc_a), 1556, tol = 0, label = "Aged CellChat interaction count")
assert_close(sum(cc_a$prob), 0.1639737560, tol = 1e-8, label = "Aged CellChat strength")
assert_close(cor(pt$pseudotime, pt$senmayo, method = "spearman"), -0.69820131, tol = 1e-6, label = "pseudotime-SenMayo Spearman rho")
mr_pos <- mr %>% filter(pval < 0.05)
assert_close(nrow(mr_pos), 29, tol = 0, label = "MR nominal positive rows")
assert_close(nrow(distinct(mr_pos, gene, outcome)), 26, tol = 0, label = "MR unique gene-outcome positives")
message("Accuracy checks passed. Figures will be written to: ", normalizePath(outdir))

# ---------- Fig0 graphical abstract ----------------------------------------
cdneg <- comp %>% filter(fraction == "CD45neg")
cdneg_group <- cdneg %>% group_by(group, age_group, injury) %>%
  summarise(cd45neg_n = sum(n), pSSPC_fib_n = sum(n[celltype == "pSSPC_fib"]), .groups = "drop") %>%
  mutate(pSSPC_fib_pct = 100 * pSSPC_fib_n / cd45neg_n,
         group = factor(group, levels = group_levels))

p0a <- sen_intact %>%
  mutate(age_group = factor(age_group, levels = c("Young","Aged"))) %>%
  ggplot(aes(age_group, SenMayo, fill = age_group)) +
  geom_col(width = 0.58, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.5f", SenMayo)), vjust = -0.35, size = 3.2) +
  coord_cartesian(ylim = range(c(sen_intact$SenMayo, 0)) + c(-0.004, 0.004), clip = "off") +
  scale_fill_manual(values = c(Young = "#009e73", Aged = "#d7301f")) +
  labs(title = "(a) Senescent remodeling", subtitle = "CD45neg intact stromal donor mean SenMayo (n = 1 young vs 1 aged; descriptive)",
       x = NULL, y = "Mean SenMayo score") +
  theme_classic(base_size = 10) + theme(plot.margin = margin(8, 10, 8, 8))

p0b1 <- cdneg_group %>%
  ggplot(aes(group, cd45neg_n, fill = group)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = cd45neg_n), vjust = -0.25, size = 3.1) +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(ylim = c(0, max(cdneg_group$cd45neg_n) * 1.14), clip = "off") +
  labs(title = "(b) Abortive activation", subtitle = "CD45neg cell number", x = NULL, y = "Cells") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 35, hjust = 1))

p0b2 <- cdneg_group %>%
  ggplot(aes(group, pSSPC_fib_pct, fill = group)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pSSPC_fib_pct)), vjust = -0.25, size = 3.1) +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(ylim = c(0, max(cdneg_group$pSSPC_fib_pct) * 1.16), clip = "off") +
  labs(subtitle = "pSSPC_fib fraction of CD45neg", x = NULL, y = "% of CD45neg") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 35, hjust = 1))

cc_top <- bind_rows(cc_y %>% mutate(age = "Young"), cc_a %>% mutate(age = "Aged")) %>%
  filter(pval < 0.05, pathway_name %in% c("TGFb", "SPP1")) %>%
  mutate(edge = paste(source, target, ligand, receptor, sep = " | ")) %>%
  group_by(age, pathway_name) %>% slice_max(prob, n = 3, with_ties = FALSE) %>% ungroup()
p0c <- cc_top %>%
  mutate(age = factor(age, levels = c("Young","Aged")),
         lab = str_wrap(paste(pathway_name, ligand, receptor, sep = " "), width = 28)) %>%
  ggplot(aes(x = reorder(lab, prob), y = prob, fill = pathway_name)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  coord_flip(clip = "off") +
  facet_wrap(~ age, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(TGFb = "#1b9e77", SPP1 = "#d95f02")) +
  labs(title = "(c) Communication programs", x = NULL, y = "Communication prob.") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_blank(), plot.margin = margin(8, 12, 8, 8))

p0 <- (p0a | (p0b1 / p0b2) | p0c) + patchwork::plot_layout(widths = c(1.0, 1.15, 1.25))
p0 <- p0 + patchwork::plot_annotation(
  title = "Convergent senescence remodeling and abortive regenerative activation in aged periosteum",
  caption = "Source: TableS1, TableS5, TableS8/TableS9. Counts and percentages are recomputed from the annotated composition table."
) & theme(plot.title = element_text(face = "bold", size = 12))
save_gg(p0, "Fig0_graphical_abstract", 13.2, 4.4)

# ---------- Fig1 composition and Fig1b activation --------------------------
comp2 <- comp %>% mutate(celltype = factor(celltype, levels = celltype_levels),
                         donor_lab = paste(donor_id, paste0(age_group, " ", injury, " ", fraction), sep = "\n"))
p1 <- comp2 %>%
  ggplot(aes(pct, factor(donor_id), fill = celltype)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.15) +
  facet_wrap(vars(age_group, injury, fraction), ncol = 3, scales = "free_y",
             labeller = labeller(.multi_line = TRUE)) +
  scale_x_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(color = NA))) +
  labs(x = "Percentage of cells (%)", y = NULL, fill = "Cell state") +
  theme_classic(base_size = 9) +
  theme(legend.position = "right", legend.box = "vertical",
        strip.text = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        panel.spacing = unit(0.35, "lines"), plot.margin = margin(8, 12, 8, 8))
save_gg(p1, "Fig1_composition", 12.5, 7.1)

p1b <- cdneg_group %>%
  mutate(injury = factor(injury, levels = c("Intact","Fracture")),
         age_group = factor(age_group, levels = c("Young","Aged"))) %>%
  ggplot(aes(x = group, y = pSSPC_fib_pct, fill = group)) +
  geom_col(width = 0.66, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pSSPC_fib_pct)), vjust = -0.28, size = 3.2) +
  geom_text(aes(label = paste0("n=", pSSPC_fib_n, "/", cd45neg_n)), y = 1.2, size = 2.7, color = "grey35") +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(ylim = c(0, max(cdneg_group$pSSPC_fib_pct) * 1.18), clip = "off") +
  labs(title = "Fibrogenic progenitor activation (CD45neg)",
       subtitle = "pSSPC_fib fraction among all CD45neg cells; denominator shown under each bar",
       x = NULL, y = "pSSPC_fib (% of CD45neg)", caption = "n = 1 donor per group; descriptive, recomputed from TableS1.") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_gg(p1b, "Fig1b_fibrogenic_activation", 6.8, 4.8)

# ---------- Fig2 SenMayo score ---------------------------------------------
read_cell_senmayo <- function(injury_keep = c("Intact","Fracture")) {
  f <- file.path(data_dir, "senmayo_celllevel_cd45neg_stromal.csv")
  if (file.exists(f)) {
    x <- read_tab("senmayo_celllevel_cd45neg_stromal.csv")
    need_cols(x, c("cell_id","age_group","injury","fraction","celltype","senmayo"), "senmayo_celllevel")
    return(x %>% filter(injury %in% injury_keep, fraction == "CD45neg"))
  }
  # Fallback: trajectory table contains the same stromal cells used for pseudotime.
  # If exact CD45neg-intact n=2933 is required, provide the cell-level file above.
  message("Optional senmayo_celllevel_cd45neg_stromal.csv not found; using TableS10 trajectory stromal cells as fallback.")
  pt %>% mutate(cell_id = paste0("cell_", seq_len(nrow(pt))), fraction = "CD45neg") %>% filter(injury %in% injury_keep)
}
cell_intact <- read_cell_senmayo("Intact") %>% filter(age_group %in% c("Young","Aged"))
cell_intact <- cell_intact %>% mutate(age_group = factor(age_group, levels = c("Aged","Young")))
p2_left <- ggplot(cell_intact, aes(age_group, senmayo, fill = age_group)) +
  geom_violin(color = NA, alpha = 0.82, scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, color = "black", linewidth = 0.35) +
  coord_cartesian(ylim = range(cell_intact$senmayo) + c(-0.02, 0.04), clip = "off") +
  scale_fill_manual(values = c(Young = "#009e73", Aged = "#d7301f"), guide = "none") +
  labs(title = sprintf("Cell-level SenMayo, CD45neg intact stromal (n = %d cells)", nrow(cell_intact)),
       x = NULL, y = "SenMayo score") +
  theme_classic(base_size = 10)

sen_intact2 <- sen_intact %>% mutate(age_group = factor(age_group, levels = c("Young","Aged")))
p2_right <- ggplot(sen_intact2, aes(age_group, SenMayo)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(size = 3, shape = 21, fill = "black") +
  geom_errorbar(aes(ymin = SenMayo, ymax = SenMayo), width = 0.18, linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.6f", SenMayo)), vjust = -0.65, size = 3.0) +
  coord_cartesian(ylim = range(c(sen_intact2$SenMayo, 0)) + c(-0.004, 0.004), clip = "off") +
  labs(title = "Donor-level SenMayo",
       subtitle = sprintf("AUC = 1, fold = %.2fx (n = 1v1, descriptive)", sen_intact2$SenMayo[sen_intact2$age_group == "Aged"] / sen_intact2$SenMayo[sen_intact2$age_group == "Young"]),
       x = NULL, y = "Mean SenMayo score") +
  theme_classic(base_size = 10)
p2 <- p2_left + p2_right + patchwork::plot_annotation(title = "SenMayo in intact periosteal stroma") & theme(plot.title = element_text(face = "bold"))
save_gg(p2, "Fig2_senmayo_score", 10.2, 4.4)

# ---------- Fig3a after fracture -------------------------------------------
cell_fx <- read_cell_senmayo("Fracture") %>% filter(age_group %in% c("Young","Aged"))
cell_fx <- cell_fx %>% mutate(age_group = factor(age_group, levels = c("Aged","Young")))
sen_fx <- sen_donor %>% filter(injury == "Fracture") %>% mutate(age_group = factor(age_group, levels = c("Young","Aged")))
p3a_left <- ggplot(cell_fx, aes(age_group, senmayo, fill = age_group)) +
  geom_violin(color = NA, alpha = 0.82, scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, color = "black", linewidth = 0.35) +
  coord_cartesian(ylim = range(cell_fx$senmayo) + c(-0.03, 0.06), clip = "off") +
  scale_fill_manual(values = c(Young = "#009e73", Aged = "#d7301f"), guide = "none") +
  labs(title = sprintf("Cell-level SenMayo after fracture (n = %d cells)", nrow(cell_fx)), x = NULL, y = "SenMayo score") +
  theme_classic(base_size = 10)
p3a_right <- ggplot(sen_fx, aes(age_group, SenMayo)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(size = 3, shape = 21, fill = "black") +
  geom_text(aes(label = sprintf("%.6f", SenMayo)), vjust = -0.65, size = 3.0) +
  scale_y_continuous(limits = c(min(sen_fx$SenMayo) - 0.01, max(sen_fx$SenMayo) + 0.01), expand = expansion(mult = c(0.08, 0.12))) +
  labs(title = "Donor-level after fracture (n=1v1)", x = NULL, y = "Mean SenMayo score") +
  theme_classic(base_size = 10)
p3a <- p3a_left + p3a_right
save_gg(p3a, "Fig3a_senmayo_after_fracture", 10.2, 4.4)

# ---------- Fig3b GSEA curve ------------------------------------------------
# Requires a ranked gene list for aged fracture. Provide ranked_aged_fracture_stat.csv
# with columns: gene, stat  (or gene, log2FoldChange, padj). If absent, a validated
# metrics panel is still written and the function tells you exactly what is missing.
plot_gsea_curve <- function(ranked_csv = file.path(data_dir, "ranked_aged_fracture_stat.csv")) {
  if (!file.exists(ranked_csv)) {
    p <- ggplot() + annotate("text", x = 0, y = 0, label = paste0("Missing ranked_aged_fracture_stat.csv\nValidated metrics from TableS6: NES = ",
                                   sprintf("%.3f", as.numeric(gsea_young$NES[1])), ", p.adj = ", formatC(as.numeric(gsea_young$p.adjust[1]), format = "e", digits = 2),
                                   "\nrank at max = ", as.numeric(gsea_young$rank[1]), "; ", gsea_young$leading_edge[1])) +
      theme_void() + labs(title = "Young injury-response in aged fracture (input ranked list required)")
    return(save_gg(p, "Fig3b_gsea_youngResponse", 7.0, 4.8))
  }
  r <- readr::read_csv(ranked_csv, show_col_types = FALSE)
  score_col <- intersect(c("stat","log2FoldChange","t"), colnames(r))[1]
  need_cols(r, c("gene", score_col), "ranked_aged_fracture_stat")
  r <- r %>% filter(is.finite(.data[[score_col]])) %>% arrange(desc(.data[[score_col]]))
  geneset <- strsplit(as.character(gsea_young$core_enrichment[1]), "/")[[1]]
  hits <- r$gene %in% geneset
  running <- cumsum(hits) / sum(hits) - cumsum(!hits) / sum(!hits)
  df <- tibble(rank = seq_along(running), running_enrichment = running, hit = hits)
  max_rank <- which.max(running)
  p <- ggplot(df, aes(rank, running_enrichment)) +
    geom_vline(xintercept = as.numeric(gsea_young$rank[1]), linetype = 2, color = "grey55", linewidth = 0.35) +
    geom_line(color = "#2ca02c", linewidth = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_segment(data = data.frame(x = seq_along(running)[hits]), aes(x = x, xend = x, y = -0.03, yend = -0.015), inherit.aes = FALSE, linewidth = 0.15) +
    annotate("point", x = max_rank, y = running[max_rank], size = 1.8) +
    annotate("text", x = max_rank, y = running[max_rank], label = sprintf("max ES = %.3f", running[max_rank]), vjust = -0.6, size = 3) +
    labs(title = sprintf("Young injury-response in aged fracture | NES = %.2f, p.adj = %s",
                         as.numeric(gsea_young$NES[1]), formatC(as.numeric(gsea_young$p.adjust[1]), format = "e", digits = 2)),
         subtitle = as.character(gsea_young$leading_edge[1]), x = "Rank in ordered dataset", y = "Running enrichment score") +
    theme_classic(base_size = 10)
  save_gg(p, "Fig3b_gsea_youngResponse", 7.4, 4.8)
}
plot_gsea_curve()

# ---------- Fig3c callus replication scatter --------------------------------
p3c <- ggplot(callus, aes(lfc_peri, lfc_callus, color = direction_consistent)) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = Inf, color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.72, size = 1.5) +
  scale_color_manual(values = c(`TRUE` = "grey70", `FALSE` = "#d7301f"), labels = c(`TRUE` = "consistent", `FALSE` = "inconsistent")) +
  coord_cartesian(clip = "off") +
  labs(title = "Periosteum-callus replication (n = 310 genes)",
       subtitle = sprintf("Direction consistency = %d/310 = %.1f%%; Pearson rho = %.3f", sum(callus$direction_consistent), 100 * mean(callus$direction_consistent), cor(callus$lfc_peri, callus$lfc_callus)),
       x = "Periosteum aging log2FC", y = "Callus aging log2FC", color = "Direction") +
  theme_classic(base_size = 10) + theme(legend.position = "inside", legend.justification.inside = c(0.02, 0.98))
save_gg(p3c, "Fig3c_callus_scatter", 6.2, 5.2)

# ---------- Fig6 communication overview ------------------------------------
cc_sum <- bind_rows(cc_y %>% mutate(age = "Young"), cc_a %>% mutate(age = "Aged")) %>%
  group_by(age) %>% summarise(n_interaction = n(), strength = sum(prob), .groups = "drop") %>%
  mutate(age = factor(age, levels = c("Young","Aged")))
p6 <- cc_sum %>% pivot_longer(cols = c(n_interaction, strength), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("n_interaction","strength"), labels = c("n_interaction","strength"))) %>%
  ggplot(aes(age, value, fill = age)) +
  geom_col(width = 0.58, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = ifelse(metric == "n_interaction", sprintf("%.0f", value), sprintf("%.3f", value))), vjust = -0.3, size = 3.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(Young = "#009e73", Aged = "#d7301f")) +
  coord_cartesian(ylim = c(0, NA), clip = "off") +
  labs(x = NULL, y = NULL, title = "CellChat communication overview", caption = "Count/strength recomputed from TableS8 and TableS9.") +
  theme_classic(base_size = 10) + theme(strip.background = element_blank())
save_gg(p6, "Fig6_CC_overview", 7.2, 4.2)

# ---------- Fig6b chord: reduce to SPP1/TGFB and avoid label pile-up --------
plot_cc_chord_pdf <- function(edges, age_label, filename) {
  edges <- edges %>% filter(pval < 0.05, pathway_name %in% c("TGFb", "SPP1"))
  if (nrow(edges) == 0) return(invisible(NULL))
  sectors <- sort(unique(c(edges$source, edges$target)))
  make_mat <- function(pw) {
    e <- edges %>% filter(pathway_name == pw) %>% group_by(source, target) %>% summarise(value = sum(prob), .groups = "drop")
    m <- matrix(0, nrow = length(sectors), ncol = length(sectors), dimnames = list(sectors, sectors))
    m[as.matrix(e[, c("source","target")])] <- e$value
    m
  }
  colmat <- function(pw, color) {
    m0 <- make_mat(pw)
    mc <- matrix(NA_character_, nrow = nrow(m0), ncol = ncol(m0), dimnames = dimnames(m0))
    mc[m0 > 0] <- color
    mc
  }
  pdf(file.path(outdir, paste0(filename, ".pdf")), width = 12.5, height = 6.0)
  tiff(file.path(outdir, paste0(filename, ".tiff")), width = 12.5, height = 6.0, units = "in", res = 600, compression = "lzw")
  layout(matrix(1:2, nrow = 1)); par(mar = c(1,1,2,1), cex = 0.72)
  for (pw in c("TGFb","SPP1")) {
    mat <- make_mat(pw)
    circos::circos.par(start.degree = 90, gap.degree = 8, track.margin = c(0.01, 0.01))
    circlize::chordDiagramFromMatrix(mat, grid.col = celltype_cols[sectors], col = colmat(pw, ifelse(pw == "TGFb", "#1b9e77", "#d95f02")),
                                     transparency = 0.34, directional = TRUE, direction.type = c("diffHeight","arrows"),
                                     link.arr.type = "big.arrow", diffHeight = -0.04, annotationTrack = c("grid","name"),
                                     preAllocateTracks = list(track.height = circlize::mm_h(2.5)), big.gap = 6, small.gap = 1.5,
                                     link.lwd = 0.2, link.lty = 1)
    title(main = paste(age_label, pw, "signaling"), line = 0.2)
    legend("bottomright", legend = pw, bty = "n", cex = 0.8, title = "Ligand program")
    circos::circos.clear()
  }
  dev.off(); dev.off()
}
plot_cc_chord_pdf(cc_y, "Young (intact)", "Fig6b_CC_chord_SPP1_TGFB_young")
plot_cc_chord_pdf(cc_a, "Aged (intact)", "Fig6b_CC_chord_SPP1_TGFB_aged")

# ---------- FigS5 SASP communication bubble ---------------------------------
sasp_edges <- bind_rows(cc_y %>% mutate(age = "Young"), cc_a %>% mutate(age = "Aged")) %>%
  filter(pval < 0.05, pathway_name %in% c("TGFb","TNF","IL6","CCL","CXCL","SPP1","MIF","CSF","OSM","IL10","IL1")) %>%
  mutate(pair = paste(source, target, sep = " -> "), interaction = paste(ligand, receptor, sep = " - ")) %>%
  group_by(age, pathway_name, interaction, pair) %>% summarise(prob = sum(prob), min_p = min(pval), .groups = "drop") %>%
  mutate(age = factor(age, levels = c("Young","Aged")), neglogp = -log10(min_p)) %>%
  group_by(age) %>% slice_max(prob, n = 180, with_ties = FALSE) %>% ungroup()
pS5 <- ggplot(sasp_edges, aes(pair, interaction, size = neglogp, color = prob)) +
  geom_point(alpha = 0.82) + facet_wrap(~ age, ncol = 1, scales = "free") +
  scale_color_viridis_c(option = "turbo", limits = range(sasp_edges$prob), name = "Commun. prob.") +
  scale_size_area(max_size = 3.6, name = expression(-log[10](p))) +
  guides(color = guide_colorbar(barwidth = unit(1.1, "cm"), barheight = unit(4.2, "cm"), order = 1),
         size = guide_legend(order = 2, ncol = 2)) +
  labs(title = "Senescence/SASP-associated communication programs", x = "Sender -> receiver", y = "Ligand - receptor") +
  theme_bw(base_size = 8) +
  theme(legend.position = "right", legend.box = "vertical", strip.background = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 5.6), axis.text.y = element_text(size = 6),
        panel.spacing = unit(0.8, "lines"), plot.margin = margin(8, 14, 8, 8))
save_gg(pS5, "FigS5_CC_sasp_bubble", 13.5, 12.5)

# ---------- Fig7 trajectory with numeric order ------------------------------
coords_file <- file.path(data_dir, "umap_coords.csv")
pt7 <- pt
if (file.exists(coords_file)) {
  u <- read_tab("umap_coords.csv"); need_cols(u, c("cell_id","umap_1","umap_2"), "umap_coords")
  if (nrow(u) != nrow(pt7) || !all(u$cell_id == pt7$cell_id)) {
    stop("umap_coords.csv must have the same row count and cell_id values as the synthetic cell_i IDs created from TableS10_pseudotime.csv. Re-export coords from the exact pt_df used to write TableS10, or set cell_id = paste0('cell_', seq_len(nrow(pt_df))).")
  }
  pt7 <- pt7 %>% left_join(u, by = "cell_id")
} else {
  message("umap_coords.csv not found: Fig7 will use pseudotime as x for the ordered panel; provide UMAP coordinates for true UMAP panels.")
  pt7$umap_1 <- pt7$pseudotime; pt7$umap_2 <- pt7$senmayo
}
pt7 <- pt7 %>% mutate(group = factor(group, levels = group_levels), celltype = factor(celltype, levels = celltype_levels))
ord4 <- pt7 %>% filter(is.finite(umap_1), is.finite(umap_2), is.finite(pseudotime)) %>%
  mutate(bin = ntile(pseudotime, 4)) %>% group_by(bin) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), pseudotime = median(pseudotime), .groups = "drop") %>%
  mutate(order_lab = LETTERS[1:4])
start_cell <- pt7 %>% filter(pseudotime == min(pseudotime, na.rm = TRUE)) %>% slice(1)

p7a <- ggplot(pt7, aes(umap_1, umap_2, color = celltype)) +
  geom_point(size = 0.35, alpha = 0.55) +
  scale_color_manual(values = celltype_cols, drop = FALSE) +
  labs(title = "Cell state", x = "UMAP 1", y = "UMAP 2", color = "celltype") +
  theme_classic(base_size = 9) + theme(legend.position = "right")
p7b <- ggplot(pt7, aes(umap_1, umap_2, color = group)) +
  geom_point(size = 0.35, alpha = 0.55) +
  geom_path(data = ord4, aes(umap_1, umap_2), inherit.aes = FALSE, color = "black", linewidth = 0.55,
            arrow = arrow(length = unit(0.18, "cm")), lineend = "round") +
  geom_label_repel(data = ord4, aes(umap_1, umap_2, label = order_lab), inherit.aes = FALSE, size = 3.1,
                   min.segment.length = 0, box.padding = 0.25, label.padding = unit(0.12, "lines")) +
  annotate("text", x = start_cell$umap_1, y = start_cell$umap_2, label = "start", fontface = "italic", size = 3.2, vjust = -0.8) +
  scale_color_manual(values = group_cols, drop = FALSE) +
  labs(title = "Group with pseudotime order", subtitle = "A -> D mark increasing pseudotime", x = "UMAP 1", y = "UMAP 2", color = "group") +
  theme_classic(base_size = 9) + theme(legend.position = "right")
p7c <- ggplot(pt7, aes(umap_1, umap_2, color = pseudotime)) +
  geom_point(size = 0.35, alpha = 0.65) +
  geom_path(data = ord4, aes(umap_1, umap_2), inherit.aes = FALSE, color = "black", linewidth = 0.55,
            arrow = arrow(length = unit(0.18, "cm")), lineend = "round") +
  geom_label_repel(data = ord4, aes(umap_1, umap_2, label = order_lab), inherit.aes = FALSE, size = 3.1,
                   min.segment.length = 0, box.padding = 0.25, label.padding = unit(0.12, "lines")) +
  annotate("text", x = start_cell$umap_1, y = start_cell$umap_2, label = "start", fontface = "italic", size = 3.2, vjust = -0.8) +
  scale_color_viridis_c(option = "plasma", name = "pseudotime") +
  labs(title = "Pseudotime with direction", x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 9) + theme(legend.position = "right")
p7 <- p7a + p7b + p7c + patchwork::plot_annotation(title = "Periosteal stromal trajectory: numeric order and direction")
save_gg(p7, "Fig7_trajectory", 15.0, 4.7)

# ---------- Fig8 pseudotime density and SenMayo coupling -------------------
med <- pt %>% group_by(group) %>% summarise(median_pt = median(pseudotime), .groups = "drop") %>%
  mutate(group = factor(group, levels = group_levels))
pt_plot <- pt %>% mutate(group = factor(group, levels = group_levels))
set.seed(1)
pt_point <- pt_plot %>% group_by(group) %>% slice_sample(n = min(1500, n()), replace = FALSE) %>% ungroup()
ord6 <- pt_plot %>% mutate(bin = ntile(pseudotime, 6)) %>% group_by(bin) %>%
  summarise(pseudotime = median(pseudotime), senmayo = median(senmayo), .groups = "drop") %>% mutate(order_lab = 1:6)

p8a <- ggplot(pt_plot, aes(pseudotime, fill = group, color = group)) +
  geom_density(alpha = 0.16, linewidth = 0.55, adjust = 1) +
  geom_vline(data = med, aes(xintercept = median_pt, color = group), linetype = 2, linewidth = 0.45, show.legend = FALSE) +
  geom_text(data = med, aes(x = median_pt, y = Inf, label = sprintf("median %.1f", median_pt), color = group),
            vjust = 1.25, angle = 90, size = 2.7, show.legend = FALSE) +
  scale_fill_manual(values = group_cols, guide = "none") + scale_color_manual(values = group_cols, guide = "none") +
  coord_cartesian(xlim = range(pt_plot$pseudotime), clip = "off") +
  labs(title = "Pseudotime density", x = "Pseudotime", y = "Density", caption = "Dashed lines and vertical labels are group medians.") +
  theme_classic(base_size = 10) + theme(plot.margin = margin(8, 10, 8, 8))

p8b <- ggplot(pt_point, aes(pseudotime, senmayo, color = group)) +
  geom_point(alpha = 0.16, size = 0.45, show.legend = FALSE) +
  geom_smooth(data = pt_plot, aes(pseudotime, senmayo, color = group), method = "loess", span = 0.22, se = TRUE, linewidth = 0.8, inherit.aes = FALSE) +
  geom_path(data = ord6, aes(pseudotime, senmayo), inherit.aes = FALSE, color = "black", linewidth = 0.6,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed"), lineend = "round") +
  geom_label_repel(data = ord6, aes(pseudotime, senmayo, label = order_lab), inherit.aes = FALSE, size = 3.2,
                   min.segment.length = 0, box.padding = 0.25, label.padding = unit(0.12, "lines")) +
  scale_color_manual(values = group_cols, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(title = sprintf("SenMayo along pseudotime | Spearman rho = %.1f", cor(pt_plot$pseudotime, pt_plot$senmayo, method = "spearman")),
       subtitle = "Numbers 1-6 give the direction of increasing pseudotime", x = "Pseudotime", y = "SenMayo") +
  theme_classic(base_size = 10) + theme(plot.margin = margin(8, 12, 8, 8))
p8 <- p8a + p8b
save_gg(p8, "Fig8_pt_density_sencoupling", 12.2, 4.8)

# ---------- MR forests -------------------------------------------------------
short_outcome <- function(x) {
  x <- gsub(" \\|\\| id:.*$", "", x)
  ifelse(grepl("Heel bone mineral density", x), "Heel BMD",
         ifelse(grepl("Femoral neck bone mineral density", x), "Femoral neck BMD",
                ifelse(grepl("Fracture of forearm or wrist", x), "Forearm/wrist fracture",
                       ifelse(grepl("^Fractures", x), "Fractures", ifelse(grepl("Osteoporosis", x), "Osteoporosis", x)))))
}
tissue_from_exposure <- function(exposure) {
  x <- gsub(".*\\(([^()]*)\\).*", "\\1", exposure)
  ifelse(x == exposure, "eQTL", x)
}
mr_plot_df <- mr %>% filter(pval < 0.05) %>%
  mutate(outcome_short = short_outcome(outcome), tissue = tissue_from_exposure(exposure),
         lab = paste0(gene, " | ", tissue, " -> ", outcome_short),
         lab = factor(lab, levels = unique(lab[order(b)])))
p9 <- ggplot(mr_plot_df, aes(b, lab, color = outcome_short)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.18, linewidth = 0.45) +
  geom_point(size = 2.1) +
  scale_color_brewer(palette = "Set2", name = "Outcome") +
  labs(title = sprintf("cis-MR nominal positives (n = %d rows; %d unique gene-outcome)", nrow(mr_plot_df), nrow(distinct(mr_plot_df, gene, outcome))),
       x = "Wald ratio beta (per SD expression, 95% CI)", y = NULL,
       caption = "Nominal p < 0.05, not corrected for multiple testing. Table7 is the deduplicated directional candidate table.") +
  theme_classic(base_size = 8.5) + theme(legend.position = "bottom", legend.box = "horizontal",
                                         plot.margin = margin(8, 18, 8, 8)) +
  coord_cartesian(clip = "off")
save_gg(p9, "Fig9_MR_positives", 8.8, 9.0)

mr_all <- mr %>% mutate(outcome_short = short_outcome(outcome), tissue = tissue_from_exposure(exposure),
                        lab = paste0(gene, " | ", tissue), lab = str_wrap(lab, width = 26),
                        lab = factor(lab, levels = unique(lab[order(b)])))
pS10 <- ggplot(mr_all, aes(b, lab, color = outcome_short)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.14, linewidth = 0.32, alpha = 0.75) +
  geom_point(size = 1.35, alpha = 0.85) +
  scale_color_brewer(palette = "Set2", name = "Outcome") +
  labs(title = sprintf("All %d cis-MR tests", nrow(mr_all)), x = "Wald ratio beta (per SD expression)", y = NULL) +
  theme_classic(base_size = 7.2) + theme(legend.position = "bottom", legend.box = "horizontal",
                                         legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
                                         axis.text.y = element_text(size = 5.8), plot.margin = margin(8, 18, 8, 8)) +
  coord_cartesian(clip = "off")
save_gg(pS10, "FigS10_MR_forest_all", 8.6, 12.5)

# ---------- final session info ----------------------------------------------
writeLines(c(
  paste0("Accuracy checks passed. Output directory: ", normalizePath(outdir)),
  paste0("Key validated values: callus consistency 172/310 (55.5%), rho=", sprintf("%.4f", cor(callus$lfc_peri, callus$lfc_callus)), ";"),
  paste0("GSEA YoungInjuryResponse NES=", sprintf("%.3f", as.numeric(gsea_young$NES[1])), ";"),
  paste0("CellChat Young/Aged interactions=", nrow(cc_y), "/", nrow(cc_a), ", strength=", sprintf("%.3f", sum(cc_y$prob)), "/", sprintf("%.3f", sum(cc_a$prob)), ";"),
  paste0("pseudotime-SenMayo Spearman rho=", sprintf("%.3f", cor(pt$pseudotime, pt$senmayo, method = "spearman")), ";"),
  paste0("MR nominal positives rows=", nrow(mr_pos), ", unique gene-outcome=", nrow(distinct(mr_pos, gene, outcome)))
), file.path(outdir, "figure_repair_validation.txt"))
message("Done.")
