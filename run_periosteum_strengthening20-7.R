#!/usr/bin/env Rscript
## ============================================================
## run_periosteum_strengthening.R v1.1
## 配套 run_periosteum_aging19.R（v3.8）——在既有输出之上系统补强三大短板：
##   短板① n=1文库无p值     → §1 bootstrap伪文库方向敏感性 + §2 置换零分布/精确检验
##   短板② SenMayo GSEA阴    → §4 多衰老签名concordance矩阵(binomial) + §5 通路级GSEA
##   短板③ 骨痂rho弱/coloc阴 → §6 跨队列meta-Z核心签名 + 尾部象限富集 + §11 coloc.susie
##   加分项                  → §3 miloR差异丰度 / §7 骨痂时间序列 / §9 多轨迹交叉验证 /
##                             §10 CellChat统计+配体活性
##
## v1.1 变更（相对 v1.0 实测反馈）：
##  A) §4: AUCell_calcAUC(list(g) → setNames(list(g), nm))，修复 "geneSets should be a named list"；
##         UCell 包缺失时自动改用 AddModuleScore 计算方向表（不再硬依赖 UCell）；
##  B) §7: 时间点解析增加 "day 7"（数字在后）正则；
##  C) §8: 新增 DO_CLOCK 开关（默认 FALSE）：GSE297256 每年龄仅 1 文库，训练 n=3 不成立，
##         时钟结果不应用于论文；确需运行则显式置 TRUE；
##  D) §2(d): SenMayo 改用中位数口径（module score 右偏，均值倍数会误导），
##         同时报告 mean fold / median fold / 中位数差置换 p；
##  E) §10: netVisual_diffInteraction 报 "replacement has length zero" → 改为手动
##         count/weight 差值热图；配体 fold 计算改为显式拆分 donor/celltype + 均值比 + 空值守卫；
##  F) §1: bootstrap 措辞收敛为 "对随机细胞抽样的敏感性"（同一 donor 对半拆非独立重复，
##         不可写成 stability validation）；
## v1.2 变更（图形修正）：
##  G) Fig14: 图例顺序 Young→Aged（breaks 显式指定）；
##  H) FigS20: 顶部 "p=" 文字被裁剪 → y 轴 expand + clip="off" + 增大上边距与图高；
##  I) FigS26: rankNet 纵坐标拥挤 → 图高 6→11、标签字号减小、边距增大；
##  J) FigS27: susie 全 NA 时补 "n.e." 文本标注、图例/坐标轴防裁剪、动态图高。
## v1.3 变更（实测反馈修复）：
##  K) §4: dir_tab 先按 signature×age_group 聚合再 pivot（v1.2 中 donor_id 参与宽表导致
##         每 donor 一行、Young/Aged 列全 NA → "0/14" 为 artifact，真实方向需重跑）；
##  L) §3: miloR buildGraph d=30 作用于 2 维 UMAP 报 "下标出界" → 改用 PCA reduced.dim, d=30；
##  M) §5: msigdbr 10.x 参数改 collection/subcollection + 504 重试 4 次；
##  N) §10: 配体活性图过滤 rec_pct<5% 的低置信受体（如 Ppbp→Cxcr2 0.57%）；
##  O) §11: runsusie 失败时把原因写入日志（解释 susie 全 NA）。
## v1.4 变更（图形修正）：
##  P) Fig11: 整图改用 ComplexHeatmap + 右侧行注释 anno_barplot(Aged−Young 有向差值),
##         替代 wrap_plots(grid.grabExpr) 拼贴（此前内容溢出/图例被截）;
##         Age 图例顺序显式 Young 在上、Aged 在下 (annotation_legend_param at=);
##  Q) FigS26: rankNet 自动过滤两组信息流均≈0/NA 的通路（如 CADM/CD226 空白行, 属并集
##         通路的正常现象）; 图高随保留通路数动态调整。
## v1.5 变更：
##  R) §5: msig_get 遍历 6 种参数组合（species="mouse"/db_species="MM" ×
##         collection/subcollection × 新旧 API），兼容各 msigdbr 版本。
## v1.6 变更：
##  S) §3: plotNhoodGraphDA 前补 buildNhoodGraph()（修复 "neighbourhood graph is missing"）;
##  T) §11: runsusie 需要 in-sample LD → 尝试经 ieugwasr::ld_matrix(pop="EUR") 获取,
##          失败维持 susie=NA 并如实报告; FigS27 图注改为 "LD unavailable"。
## v1.7 变更（图形修正）：
##  V) FigS22: age_group 显式 factor(Young,Aged) → Young 左 Aged 右；
##  W) FigS23b: 图幅 7×4.5 → 10×6.5, 通路名去 REACTOME_/HALLMARK_ 前缀并换行,
##             x 轴加 expand 缓解柱压缩；
##  X) FigS21: 火山图按 logFC 方向着色; 邻域图副标题注明 "0 显著 → 全灰属预期 (n=1v1)"。
## 运行: Rscript run_periosteum_strengthening.R（主脚本 v3.8 之后、同一工作目录）
## ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(data.table); library(Matrix)
  library(patchwork); library(ComplexHeatmap); library(clusterProfiler)
  library(org.Mm.eg.db); library(cowplot)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; count <- dplyr::count

dir.create("output/figures_strength", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
log_file <- sprintf("output/logs/strengthening_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " "))
  cat(msg, "\n"); write(msg, log_file, append = TRUE)
}
run_step <- function(nm, expr) {
  log_msg("---- STEP: ", nm)
  tryCatch(eval(substitute(expr), envir = parent.frame()),
           error = function(e) { log_msg("!! ERROR in ", nm, ": ", e$message); NULL })
}
save_fig <- function(name, width, height, expr, dpi = 600,
                     dir = "output/figures_strength") {
  e <- substitute(expr); envir <- parent.frame()
  pdf(file.path(dir, paste0(name, ".pdf")), width = width, height = height)
  print(eval(e, envir)); dev.off()
  tiff(file.path(dir, paste0(name, ".tiff")), width = width, height = height,
       units = "in", res = dpi, compression = "lzw", type = "cairo")
  print(eval(e, envir)); dev.off()
  log_msg("图已生成: ", name)
}
spearman_rho <- function(x, y) {
  r <- suppressWarnings(tryCatch(cor.test(x, y, method = "spearman"),
                                 error = function(e) NULL))
  if (is.null(r)) NA_real_ else unname(r$estimate)
}
safe_csv <- function(path) if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL

## ---------------- 全局配置 ----------------
PERI_RDS    <- "output/peri_seurat.rds"
N_BOOT      <- 100      # §1 bootstrap 迭代数
N_PERM      <- 1000     # §2 细胞级置换
N_PERM_BIG  <- 10000    # §2 表格级置换
DO_MILO     <- TRUE
DO_TRADESEQ <- FALSE    # 计算量大, 确认后改 TRUE
DO_AUCELL   <- TRUE
DO_CLOCK    <- FALSE    # v1.1: 训练集 n=3, 时钟不用于论文；确需运行改 TRUE
STROMAL     <- c("pSSPC_q","pSSPC_fib","pSSPC_osteo","Endo","Cycling","Chondro")
SENDERS     <- c("Neutroph","Macroph","Endo","pSSPC_fib","Chondro")
AGE_MONTHS  <- c("Young" = 3.5, "Aged" = 22)

## ============================================================
## §1 细胞级 bootstrap 伪文库 + 方向敏感性（针对 n=1 文库）
## ============================================================
run_step("bootstrap_stability", {
  if (!file.exists(PERI_RDS)) {
    log_msg("缺 ", PERI_RDS, " → 跳过"); return(NULL) }
  peri <- readRDS(PERI_RDS)
  peri_intact <- subset(peri, subset = injury == "Intact" & fraction == "CD45neg" &
                          celltype %in% c("pSSPC_q","pSSPC_fib","pSSPC_osteo","Endo","Chondro"))
  rm(peri); gc(reset = TRUE)
  counts <- GetAssayData(peri_intact, assay = "RNA", layer = "counts")
  meta <- peri_intact@meta.data

  agg_once <- function(counts, meta, seed) {
    set.seed(seed)
    half <- unsplit(lapply(split(seq_len(nrow(meta)), meta$donor_id),
                           function(i) sample(rep(1:2, length.out = length(i)))),
                    meta$donor_id)
    plib <- paste(meta$donor_id, half, sep = "_")
    plibs <- unique(plib)
    M <- do.call(cbind, lapply(plibs, function(p)
      Matrix::rowSums(counts[, plib == p, drop = FALSE])))
    colnames(M) <- plibs
    lcpm <- log2(sweep(M, 2, colSums(M), "/") * 1e6 + 1)
    age <- meta$age_group[match(plibs, plib)]
    setNames(rowMeans(lcpm[, age == "Aged", drop = FALSE]) -
             rowMeans(lcpm[, age == "Young", drop = FALSE]), rownames(M))
  }
  log_msg("bootstrap 伪文库 DE: ", N_BOOT, " 轮迭代 ...")
  lfc_mat <- vapply(seq_len(N_BOOT), function(i) agg_once(counts, meta, 1000 + i),
                    numeric(nrow(counts)))
  rownames(lfc_mat) <- rownames(counts)
  keep <- rowMeans(lfc_mat != 0) >= 0.5
  lfc_mat <- lfc_mat[keep, , drop = FALSE]

  lfc_med  <- apply(lfc_mat, 1, median)
  dir_prob <- rowMeans(sign(lfc_mat) == sign(lfc_med))
  stab <- tibble(gene = rownames(lfc_mat), lfc_median = lfc_med,
                 dir_prob = dir_prob, stable = dir_prob >= 0.8) %>%
    arrange(desc(abs(lfc_median)))
  write.csv(stab, "output/tables/TableS12_stability_bootstrap.csv", row.names = FALSE)

  deg1 <- safe_csv("output/tables/Table1_aging_DEGs.csv")
  stab_pct <- NA_real_
  if (!is.null(deg1)) {
    ov <- stab %>% filter(gene %in% deg1$gene)
    stab_pct <- round(mean(ov$stable) * 100, 1)
    log_msg("方向敏感性: Table1 DEGs 中 ", stab_pct, "% 在 ", N_BOOT,
            " 轮随机细胞抽样中方向一致")
    log_msg("注意: 同一 donor 细胞对半拆分不构成独立重复, 本分析仅说明对抽样不敏感,",
            " 不可表述为 stability validation")
  }

  top_show <- stab %>% filter(abs(lfc_median) > 0.15) %>%
    arrange(dir_prob) %>% slice_head(n = 4) %>%
    bind_rows(stab %>% filter(abs(lfc_median) > 0.15, stable) %>%
                arrange(desc(abs(lfc_median))) %>% slice_head(n = 4)) %>%
    distinct(gene) %>% pull(gene)
  plot_df <- as.data.frame(t(lfc_mat[top_show, , drop = FALSE])) %>%
    pivot_longer(everything(), names_to = "gene", values_to = "lfc_iter") %>%
    left_join(stab %>% select(gene, dir_prob), by = "gene") %>%
    mutate(gene = sprintf("%s (P_same=%.2f)", gene, dir_prob))
  p_dist <- ggplot(plot_df, aes(lfc_iter, reorder(gene, dir_prob))) +
    ggridges::geom_density_ridges(scale = 1.1, rel_min_height = 0.01,
                                  fill = "#4DBBD5", alpha = .5) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    labs(x = "log2FC per iteration (Aged - Young, pseudo-bulk 2v2)",
         y = NULL, title = "Sensitivity to random cell subsampling",
         subtitle = paste0("each donor split into 2 pseudo-libraries, n = ", N_BOOT,
                           " iterations (non-independent; sensitivity analysis)"))
  p_bar <- ggplot(stab %>% filter(abs(lfc_median) > 0.1), aes(dir_prob)) +
    geom_histogram(binwidth = 0.05, fill = "#00A087", color = "white") +
    geom_vline(xintercept = 0.8, linetype = 2, color = "#E64B35") +
    labs(x = paste0("P(same direction) over ", N_BOOT, " iterations"),
         y = "Genes", title = "Direction-probability distribution",
         caption = "pseudo-libraries derived from the same donor are not independent replicates")
  save_fig("Fig10_bootstrap_stability", 12.5, 5.2,
           print(wrap_plots(p_dist, p_bar, ncol = 2)))
  rm(lfc_mat, counts, peri_intact); gc(reset = TRUE)
})

## ============================================================
## §2 置换零分布 + 精确检验（给描述性数字配实证 p）
## ============================================================
run_step("permutation_nulls", {
  emp <- list()

  cv <- safe_csv("output/tables/Table4_callus_validation.csv")
  if (!is.null(cv) && nrow(cv) > 50) {
    rho_obs <- spearman_rho(cv$lfc_peri, cv$lfc_callus)
    perm <- replicate(N_PERM_BIG, spearman_rho(cv$lfc_peri, sample(cv$lfc_callus)))
    perm <- perm[is.finite(perm)]
    emp$callus_rho <- tibble(
      metric = "callus_spearman_rho_global", observed = round(rho_obs, 3),
      null_median = round(median(perm), 3),
      empirical_p = signif(mean(abs(perm) >= abs(rho_obs)), 3),
      n_perm = length(perm),
      note = "label permutation; global correlation diluted -> interpret with tail-quadrant test (TableS17c)")
    bt <- binom.test(sum(cv$direction_consistent, na.rm = TRUE),
                     sum(!is.na(cv$direction_consistent)), 0.5)
    emp$callus_consistency <- tibble(
      metric = "callus_direction_consistency", observed = round(mean(cv$direction_consistent, na.rm=TRUE), 3),
      null_median = 0.5, empirical_p = signif(bt$p.value, 3), n_perm = NA_integer_,
      note = paste0("exact binomial test, n=", sum(!is.na(cv$direction_consistent))))
  }

  cv2 <- safe_csv("output/tables/Table4d_atlas_validation.csv")
  if (!is.null(cv2) && nrow(cv2) > 50) {
    rho_obs2 <- spearman_rho(cv2$lfc_peri, cv2$slope_atlas)
    perm2 <- replicate(N_PERM_BIG, spearman_rho(cv2$lfc_peri, sample(cv2$slope_atlas)))
    perm2 <- perm2[is.finite(perm2)]
    emp$atlas_rho <- tibble(
      metric = "atlas_spearman_rho", observed = round(rho_obs2, 3),
      null_median = round(median(perm2), 3),
      empirical_p = signif(mean(abs(perm2) >= abs(rho_obs2)), 3),
      n_perm = length(perm2), note = "label permutation of atlas slope")
    bt2 <- binom.test(sum(cv2$direction_consistent, na.rm = TRUE),
                      sum(!is.na(cv2$direction_consistent)), 0.5)
    emp$atlas_consistency <- tibble(
      metric = "atlas_direction_consistency",
      observed = round(mean(cv2$direction_consistent, na.rm = TRUE), 3),
      null_median = 0.5, empirical_p = signif(bt2$p.value, 3),
      n_perm = NA_integer_, note = "exact binomial test")
  }

  pt <- safe_csv("output/tables/TableS10_pseudotime.csv")
  if (!is.null(pt)) {
    med <- function(g) median(pt$pseudotime[pt$group == g], na.rm = TRUE)
    d_obs <- (med("Young_Fracture") - med("Young_Intact")) -
             (med("Aged_Fracture")  - med("Aged_Intact"))
    grp_all <- pt$group
    perm_d <- replicate(N_PERM, {
      g2 <- sample(grp_all)
      m2 <- function(gg) median(pt$pseudotime[g2 == gg], na.rm = TRUE)
      (m2("Young_Fracture") - m2("Young_Intact")) -
      (m2("Aged_Fracture")  - m2("Aged_Intact"))
    })
    boot_d <- replicate(N_PERM, {
      i <- sample(seq_len(nrow(pt)), replace = TRUE)
      ptb <- pt[i, ]
      mb <- function(gg) median(ptb$pseudotime[ptb$group == gg], na.rm = TRUE)
      (mb("Young_Fracture") - mb("Young_Intact")) -
      (mb("Aged_Fracture")  - mb("Aged_Intact"))
    })
    emp$pt_double_diff <- tibble(
      metric = "pseudotime_progression_double_diff", observed = round(d_obs, 3),
      null_median = round(median(perm_d), 3),
      empirical_p = signif(mean(abs(perm_d) >= abs(d_obs)), 3),
      n_perm = N_PERM,
      note = paste0("bootstrap 95% CI [", round(quantile(boot_d, .025), 2), ", ",
                    round(quantile(boot_d, .975), 2), "]; cell-level, exploratory"))
    log_msg("Δpseudotime 双重差 = ", round(d_obs, 2), " | 置换 p = ",
            signif(mean(abs(perm_d) >= abs(d_obs)), 3))
    r_obs <- spearman_rho(pt$pseudotime, pt$senmayo)
    perm_r <- replicate(N_PERM, spearman_rho(pt$pseudotime[sample(seq_len(nrow(pt)))],
                                             pt$senmayo))
    perm_r <- perm_r[is.finite(perm_r)]
    emp$pt_senmayo_rho <- tibble(
      metric = "pseudotime_x_senmayo_rho", observed = round(r_obs, 3),
      null_median = round(median(perm_r), 3),
      empirical_p = signif(mean(abs(perm_r) >= abs(r_obs)), 3),
      n_perm = length(perm_r), note = "joint permutation, cell-level")
  }

  ## v1.1(d): SenMayo 用中位数口径（module score 右偏，均值比会误导）
  if (file.exists(PERI_RDS)) {
    peri <- readRDS(PERI_RDS)
    sc <- peri@meta.data %>% filter(fraction == "CD45neg", injury == "Intact")
    yv <- sc$SenMayo1[sc$age_group == "Young"]; av <- sc$SenMayo1[sc$age_group == "Aged"]
    d_obs <- median(av) - median(yv)
    allv <- c(yv, av); ny <- length(yv)
    perm_d <- replicate(N_PERM, {
      s <- sample(allv); median(s[seq_len(ny)]) - median(s[-seq_len(ny)]) })
    emp$senmayo_median_diff <- tibble(
      metric = "senmayo_median_diff_aged_minus_young", observed = round(d_obs, 4),
      null_median = round(median(perm_d), 4),
      empirical_p_greater = signif(mean(perm_d >= d_obs), 3),
      empirical_p_two = signif(mean(abs(perm_d) >= abs(d_obs)), 3),
      n_perm = N_PERM,
      note = paste0("median fold = ", round(median(av)/median(yv), 3),
                    "; mean fold = ", round(mean(av)/mean(yv), 3),
                    " (right-skewed; median preferred)"))
    rm(peri); gc(reset = TRUE)
  }

  emp_df <- bind_rows(emp)
  write.csv(emp_df, "output/tables/TableS13_empirical_pvalues.csv", row.names = FALSE)
  print(as.data.frame(emp_df))
  ## v1.2(H): expand y 轴 + clip="off" + 增大上/右边距, 防止顶部 "p=" 标签被裁剪
  save_fig("FigS20_empirical_p_forest", 9, max(4.5, 0.5 * nrow(emp_df) + 2), {
  print(ggplot(emp_df, aes(observed, reorder(metric, observed))) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
        geom_point(size = 3, color = "#3C5488") +
        geom_text(aes(label = paste0("p=", empirical_p)), vjust = -0.55, size = 3,
                  color = "#E64B35", lineheight = .9) +
        scale_y_discrete(expand = expansion(mult = c(0.18, 0.18))) +
        coord_cartesian(clip = "off") +
        labs(x = "Observed statistic", y = NULL,
             title = "Empirical null-distribution p values",
             caption = "exact binomial for consistency; permutation otherwise") +
        theme(plot.title = element_text(size = 11, face = "bold"),
              plot.margin = margin(14, 20, 6, 6)))
  })
})

## ============================================================
## §3 miloR 邻域差异丰度（组成变化的图统计升级）
## ============================================================
run_step("miloR", {
  if (!DO_MILO) return(NULL)
  if (!file.exists(PERI_RDS)) { log_msg("缺 peri RDS → 跳过 miloR"); return(NULL) }
  if (!requireNamespace("miloR", quietly = TRUE)) { log_msg("缺 miloR 包"); return(NULL) }
  peri <- readRDS(PERI_RDS)
  sce <- as.SingleCellExperiment(peri)
  SingleCellExperiment::reducedDim(sce, "UMAP") <- Embeddings(peri, "umap")
  milo <- miloR::Milo(sce)
  ## v1.3(L): UMAP 仅 2 维, d=30 会 "下标出界"; 图构建用 PCA(d=30), 出图仍用 UMAP
  milo <- miloR::buildGraph(milo, k = 20, d = 30, reduced.dim = "PCA")
  milo <- miloR::makeNhoods(milo, prop = 0.1, k = 20, d = 30, refined = TRUE)
  milo <- miloR::countCells(milo, meta.data = as.data.frame(SummarizedExperiment::colData(milo)),
                            samples = "donor_id")
  design.df <- as.data.frame(SummarizedExperiment::colData(milo)) %>%
    distinct(donor_id, age_group) %>% as.data.frame()
  rownames(design.df) <- design.df$donor_id
  milo_res <- miloR::testNhoods(milo, design = ~ age_group, design.df = design.df,
                                fdr.weighting = "k-distance")
  ## v1.6(S): plotNhoodGraphDA 需要邻域图, 否则报 "neighbourhood graph is missing"
  milo <- miloR::buildNhoodGraph(milo)
  log_msg("miloR: FDR<0.05 邻域数 = ", sum(milo_res$SpatialFDR < 0.05, na.rm = TRUE),
          " / ", nrow(milo_res))
  write.csv(milo_res, "output/tables/TableS14_milo_nhoods.csv", row.names = FALSE)
  ## v1.7(X): 0 显著邻域 → 全灰属预期; 火山图按方向着色更有信息量
  save_fig("FigS21_milo_DA", 12, 5, {
    print(wrap_plots(
      miloR::plotNhoodGraphDA(milo, milo_res, layout = "UMAP", alpha = 0.1) +
        ggtitle("Milo neighbourhoods (FDR<0.1)") +
        labs(subtitle = paste0(sum(milo_res$SpatialFDR < 0.1, na.rm = TRUE),
                               " significant — uniform grey is the expected null at n=1v1")),
      ggplot(milo_res, aes(logFC, -log10(SpatialFDR), color = logFC > 0)) +
        geom_point(size = .5, alpha = .5) +
        scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#00A087"),
                           labels = c("TRUE" = "Aged", "FALSE" = "Young")) +
        geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey30") +
        labs(x = "logFC Aged vs Young", y = expression(-log[10]~FDR),
             title = "Neighbourhood DA volcano", color = "Enriched in"),
      ncol = 2))
  })
  rm(peri, sce, milo); gc(reset = TRUE)
})

## ============================================================
## §4 多衰老签名 concordance 矩阵（针对 SenMayo GSEA 阴性）
## ============================================================
run_step("multi_signature", {
  if (!file.exists(PERI_RDS)) { log_msg("缺 peri RDS → 跳过签名矩阵"); return(NULL) }
  SEN_SIGNATURES <- list(
    SenMayo_core = c("Ccl2","Il6","Serpine1","Icam1","Timp1","Mmp10","Cxcl8","Plau",
                     "Il1b","Timp2","Igfbp7","Mmp1a","Mmp3","Tnf","Fas","Ccl20","Ccl3",
                     "Ccl5","Cxcl1","Cxcl2","Spp1","Angptl4","Gdf15","Igfbp3","Igfbp5",
                     "Igfbp6","Pappa","Stc1","Mif","Hmgb1","Hmgb2","Anxa1","Vamp3",
                     "Vamp5","Vamp7","Arhgdib","Capg","Cops5"),
    CellAge_core_up = c("Cdkn2a","Trp53","Cdkn1a","Trp21","Serpine1","Tgfb1","Igfbp3",
                        "Igfbp5","Mmp3","Il6","Ccl2","Fos","Jun","Egr1","Ddit4","Gdf15",
                        "Plau","Vegfa","Timp1","Icam1","Cxcl1","Fgf2","Pdgfb","Serpine2",
                        "Mif","Hmgb1","Btg1","Gadd45a","Sesn1","Txnip"),
    Fridman_core = c("Il6","Serpine1","Tgfb1","Mmp3","Igfbp3","Igfbp5","Plau","Vegfa",
                     "Fos","Jun","Egr1","Btg1","Gadd45a","Ccl2","Icam1","Timp1","F3",
                     "Serpine2","Mif"),
    SASP_Coppe = c("Il6","Cxcl1","Cxcl2","Ccl2","Ccl5","Mmp3","Mmp10","Serpine1","Vegfa",
                   "Tnf","Icam1","Il1b","Csf2","Csf3","Gdf15","Spp1","Il1a"),
    DecipherS = c("Cdkn2a","Cdkn1a","Trp53","Trp21","Bax","Serpine1","Tgfb1","Igfbp7",
                  "Fos","Jun","Egr1","Mxi1","Zmat3","Sertad1","Ddit4","Rprml","Gadd45a"),
    LM_aging = c("Cdkn2a","Cdkn1a","Trp53","Trp21","Gadd45a","Sesn1","Ddit4","Txnip",
                 "Btg1","Zmat3","Sod2","Cat","Gsr","Nqo1"),
    GOBP_senescence = c("Cdkn1a","Cdkn2a","Trp53","Trp21","Gadd45a","Bax","Bcl2","Mdm2",
                        "Casp3","Mapk14","Nfkb1","Tnf","Il6","Fas")
  )
  peri <- readRDS(PERI_RDS)

  ## 方法1: UCell（若已安装）
  ucell_ok <- requireNamespace("UCell", quietly = TRUE)
  if (ucell_ok) {
    sigs_ok <- lapply(SEN_SIGNATURES, function(s) intersect(s, rownames(peri)))
    sigs_ok <- sigs_ok[lengths(sigs_ok) >= 5]
    peri <- UCell::AddModuleScore_UCell(peri, features = sigs_ok, name = "_UCell")
    log_msg("UCell 签名: ", paste(names(sigs_ok), collapse = ", "))
  }
  ## 方法2: Seurat AddModuleScore
  for (nm in names(SEN_SIGNATURES)) {
    g <- intersect(SEN_SIGNATURES[[nm]], rownames(peri))
    if (length(g) >= 5)
      peri <- AddModuleScore(peri, features = list(g),
                             name = paste0("SIG_", nm, "_"))
  }
  ## 方法3: AUCell（子采样；v1.1: 命名 list 修复）
  auc_cols <- character(0)
  if (DO_AUCELL && requireNamespace("AUCell", quietly = TRUE)) {
    set.seed(7)
    cs <- colnames(peri)[sample(ncol(peri), min(20000, ncol(peri)))]
    expr_rank <- AUCell::AUCell_buildRankings(
      GetAssayData(peri, assay = "RNA", layer = "counts")[, cs])
    for (nm in names(SEN_SIGNATURES)) {
      g <- intersect(SEN_SIGNATURES[[nm]], rownames(peri))
      if (length(g) >= 5) {
        auc <- AUCell::AUCell_calcAUC(setNames(list(g), nm), expr_rank)
        peri@meta.data[cs, paste0("AUC_", nm)] <- as.numeric(AUCell::getAUC(auc)[1, ])
        auc_cols <- c(auc_cols, paste0("AUC_", nm))
      }
    }
  }
  gc(reset = TRUE)

  md <- peri@meta.data %>% filter(fraction == "CD45neg")
  sig_cols <- c(grep("_UCell$", colnames(md), value = TRUE),
                grep("^SIG_.*_1$", colnames(md), value = TRUE), auc_cols)
  donor_sig <- md %>%
    select(donor_id, age_group, injury, celltype, all_of(sig_cols)) %>%
    pivot_longer(all_of(sig_cols), names_to = "sig_col", values_to = "score") %>%
    mutate(method = case_when(str_detect(sig_col, "_UCell$") ~ "UCell",
                              str_detect(sig_col, "^SIG_")   ~ "AddModuleScore",
                              TRUE ~ "AUCell"),
           signature = sig_col %>%
             str_remove("_UCell$") %>% str_remove("^SIG_") %>% str_remove("_1$") %>%
             str_remove("^AUC_")) %>%
    group_by(donor_id, age_group, injury, method, signature) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  write.csv(donor_sig, "output/tables/TableS15_multisig_donor.csv", row.names = FALSE)

  ## v1.3(K): 先按 signature×age_group 聚合再宽表。
  ## v1.2 的写法把 donor_id 带进了 pivot_wider → 每个 donor 单独一行,
  ## Young/Aged 两列几乎全 NA → up_in_aged 恒为 NA → "0/14" 是 artifact, 必须重跑本步。
  method_main <- if (ucell_ok && "UCell" %in% donor_sig$method) "UCell" else "AddModuleScore"
  dir_tab <- donor_sig %>% filter(method == method_main, injury == "Intact") %>%
    group_by(signature, age_group) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = age_group, values_from = score) %>%
    mutate(up_in_aged = Aged > Young)
  n_up <- sum(dir_tab$up_in_aged, na.rm = TRUE); n_sig <- nrow(dir_tab)
  binom_p <- binom.test(n_up, n_sig, 0.5)$p.value
  log_msg("多签名 concordance(", method_main, "): ", n_up, "/", n_sig,
          " 个签名在老龄骨膜升高 | 二项检验 p = ", signif(binom_p, 3))

  meth_tab <- donor_sig %>% filter(injury == "Intact") %>%
    group_by(signature, method) %>%
    summarise(up = sum(score[age_group == "Aged"], na.rm = TRUE) >
                     sum(score[age_group == "Young"], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = method, values_from = up)
  meth_tab$agree <- rowSums(as.matrix(meth_tab[, -1]), na.rm = TRUE)
  write.csv(meth_tab, "output/tables/TableS15b_multisig_method_concordance.csv",
            row.names = FALSE)

  ## v1.4(P): 整图为单个 ComplexHeatmap 对象: 热图 + 右侧 anno_barplot 行注释。
  ## 行注释数据 = dir_tab 的 Aged − Young 供体分数差（即原 "Aged > Young ?" 条的
  ## 数值化, 数据来源 TableS15_multisig_donor.csv, method_main × Intact）。
  heat_df <- donor_sig %>% filter(injury == "Intact")
  mat <- heat_df %>% filter(method == method_main) %>%
    select(signature, donor_id, score) %>%
    pivot_wider(names_from = donor_id, values_from = score) %>%
    column_to_rownames("signature") %>% as.matrix()
  mat_z <- t(scale(t(mat)))
  hm_meta <- heat_df %>% filter(method == method_main) %>%
    distinct(donor_id, age_group) %>%
    arrange(match(donor_id, colnames(mat_z)))
  delta <- dir_tab$Aged - dir_tab$Young
  names(delta) <- dir_tab$signature
  ord <- rownames(mat_z)[order(delta[rownames(mat_z)])]
  mat_z <- mat_z[ord, , drop = FALSE]
  ht <- Heatmap(mat_z, name = "z-score",
                top_annotation = HeatmapAnnotation(
                  Age = hm_meta$age_group,
                  col = list(Age = c("Young" = "#00A087", "Aged" = "#E64B35")),
                  annotation_legend_param = list(Age = list(at = c("Young", "Aged")))),
                right_annotation = rowAnnotation(
                  `Aged - Young` = anno_barplot(
                    delta[rownames(mat_z)],
                    gp = gpar(fill = ifelse(delta[rownames(mat_z)] > 0, "#E64B35", "#00A087")),
                    width = unit(4.2, "cm"),
                    axis_param = list(gp = gpar(fontsize = 7)))),
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_names_gp = gpar(fontsize = 8.5),
                column_names_gp = gpar(fontsize = 8),
                column_title = paste0(method_main, " scores (CD45neg intact, donor level; ",
                                      n_up, "/", n_sig, " up in aged, binom p = ",
                                      signif(binom_p, 2), ")"))
  save_fig("Fig11_multisig_heatmap", 12.5, 5.8, { draw(ht) })
  ## v1.7(V): Young 左 Aged 右
  save_fig("FigS22_multisig_methods", 10, 4.5, {
    print(ggplot(donor_sig %>% filter(injury == "Intact") %>%
                   mutate(age_group = factor(age_group, levels = c("Young", "Aged"))),
                 aes(age_group, score, fill = age_group)) +
          geom_boxplot(outlier.size = .4) +
          facet_grid(method ~ signature, scales = "free_y") +
          scale_fill_manual(values = c("Young" = "#00A087", "Aged" = "#E64B35")) +
          theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 7)) +
          labs(x = NULL, y = "Signature score", title = "Method robustness across scoring algorithms"))
  })
  rm(peri); gc(reset = TRUE)
})

## ============================================================
## §5 通路级 GSEA（HALLMARK + REACTOME，含衰老/炎症轴）
## ============================================================
run_step("pathway_GSEA", {
  f <- "output/tables/TableS0_full_DE_stats.csv"
  if (!file.exists(f)) { log_msg("缺 ", f, " → 跳过通路GSEA"); return(NULL) }
  if (!requireNamespace("msigdbr", quietly = TRUE)) { log_msg("缺 msigdbr"); return(NULL) }
  res <- read.csv(f, stringsAsFactors = FALSE)
  gene_list <- if ("stat" %in% names(res) && all(!is.na(res$stat))) {
    res$stat
  } else sign(res$log2FoldChange) * -log10(res$pvalue + 1e-300)
  names(gene_list) <- res$gene
  gene_list <- sort(gene_list, decreasing = TRUE)

  ## v1.5(R): msigdbr 各版本 API 差异大（species/db_species × collection/category ×
  ##          subcollection/subcategory）→ 遍历 6 种组合, 成功的第一套即采用
  msig_get <- function(collection, subcollection = NULL) {
    variants <- list(
      function() msigdbr::msigdbr(species = "mouse", collection = collection,
                                  subcollection = subcollection),
      function() msigdbr::msigdbr(species = "mouse", collection = collection,
                                  subcategory  = subcollection),
      function() msigdbr::msigdbr(db_species = "MM", collection = collection,
                                  subcollection = subcollection),
      function() msigdbr::msigdbr(db_species = "MM", collection = collection,
                                  subcategory  = subcollection),
      function() msigdbr::msigdbr(species = "Mus musculus", category = collection,
                                  subcategory = subcollection),
      function() msigdbr::msigdbr(db_species = "MM", category = collection,
                                  subcategory = subcollection)
    )
    for (i in seq_along(variants)) {
      out <- tryCatch(variants[[i]](), error = function(e) {
        log_msg("  msigdbr 尝试", i, "失败: ", e$message); NULL })
      if (!is.null(out) && nrow(out) > 0) {
        log_msg("  msigdbr 成功(方案", i, ")")
        return(out)
      }
      Sys.sleep(3)
    }
    NULL
  }
  mH  <- msig_get("H")
  mC2 <- msig_get("C2", "CP:REACTOME")
  if (is.null(mH) || is.null(mC2)) {
    log_msg("msigdbr 获取失败（网络/版本）→ 本步跳过，可稍后单独重跑"); return(NULL) }
  t2g <- bind_rows(mH, mC2) %>% distinct(gs_name, gene_symbol) %>%
    rename(term = gs_name, gene = gene_symbol)
  set.seed(123)
  g <- GSEA(gene_list, TERM2GENE = t2g, pAdjustMethod = "BH",
            minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, eps = 1e-10)
  gdf <- g@result
  write.csv(gdf, "output/tables/TableS16_pathway_GSEA.csv", row.names = FALSE)
  sen_axis <- gdf %>% filter(str_detect(tolower(ID),
      "senescen|sasp|aging|ageing|p53|nfkb|inflam|tnfa|il6|jak|interferon|apopt"))
  write.csv(sen_axis, "output/tables/TableS16b_senescence_axis_GSEA.csv", row.names = FALSE)
  log_msg("通路GSEA: 显著(FDR<0.05) = ", sum(gdf$p.adjust < 0.05),
          "；衰老/炎症轴显著 = ", sum(sen_axis$p.adjust < 0.05))

  top_up <- gdf %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = 15)
  top_dn <- gdf %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = 15)
  save_fig("FigS23_pathway_GSEA", 9.5, 8, {
    print(ggplot(bind_rows(top_up, top_dn),
                 aes(NES, reorder(ID, NES), size = -log10(p.adjust),
                     color = p.adjust < 0.05)) +
          geom_point() + scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey60")) +
          geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
          labs(x = "NES", y = NULL, size = expression(-log[10]~FDR), color = "FDR<0.05",
               title = "Pathway-level GSEA of periosteal aging",
               subtitle = "HALLMARK + REACTOME (full stats: TableS16)") +
          theme(legend.position = "top", axis.text.y = element_text(size = 7)))
  })
  if (nrow(sen_axis) > 0) {
    ## v1.7(W): 去前缀+换行缩短标签, 大图幅, x 轴 expand 缓解压缩
    sa_plot <- sen_axis %>% arrange(p.adjust) %>% slice_head(n = 12) %>%
      mutate(ID = stringr::str_wrap(stringr::str_remove(ID, "^REACTOME_|^HALLMARK_"),
                                    width = 38))
    save_fig("FigS23b_senescence_axis", 10, 6.5, {
      print(ggplot(sa_plot, aes(NES, reorder(ID, NES), fill = p.adjust < 0.05)) +
            geom_col() +
            scale_x_continuous(expand = expansion(mult = c(0.08, 0.18))) +
            scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey60")) +
            labs(x = "NES", y = NULL, fill = "FDR<0.05",
                 title = "Senescence / p53-apoptosis / inflammatory axis") +
            theme(axis.text.y = element_text(size = 7.5),
                  plot.margin = margin(6, 14, 6, 6)))
    })
  }
})

## ============================================================
## §6 跨队列 meta-Z 核心签名 + 骨痂尾部象限富集（针对 rho≈0.05）
## ============================================================
run_step("cross_cohort_meta", {
  mk_z <- function(df, lfc_col, p_col) {
    df %>% filter(!is.na(.data[[p_col]]), .data[[p_col]] > 0,
                  !is.na(.data[[lfc_col]])) %>%
      mutate(z = sign(.data[[lfc_col]]) * abs(qnorm(.data[[p_col]] / 2))) %>%
          select(gene, z, lfc = all_of(lfc_col))
  }
  zl <- list()
  f0 <- "output/tables/TableS0_full_DE_stats.csv"
  if (file.exists(f0)) zl$periosteum <- mk_z(read.csv(f0), "log2FoldChange", "pvalue")
  f4b <- "output/tables/Table4b_callus_aging_DEGs.csv"
  if (file.exists(f4b)) zl$callus <- mk_z(read.csv(f4b), "log2FoldChange", "pvalue")
  f4c <- "output/tables/Table4c_atlas_age_slope.csv"
  if (file.exists(f4c)) zl$atlas <- mk_z(read.csv(f4c), "log2FoldChange", "pvalue")
  zl <- compact(zl)
  if (length(zl) < 2) { log_msg("可用队列<2 → 跳过 meta-Z"); return(NULL) }
  log_msg("meta-Z 队列: ", paste(names(zl), collapse = " + "))

  zs <- bind_rows(zl, .id = "cohort")
  meta <- zs %>% group_by(gene) %>%
    summarise(z_meta = sum(z) / sqrt(n()), n_cohort = n(),
              dir_consistent = abs(sum(sign(z))) == n(), .groups = "drop") %>%
    mutate(p_meta = 2 * pnorm(-abs(z_meta)),
           fdr = p.adjust(p_meta, "BH")) %>%
    arrange(p_meta)
  write.csv(meta, "output/tables/TableS17_metaZ.csv", row.names = FALSE)
  core <- meta %>% filter(dir_consistent, n_cohort >= 3, fdr < 0.05)
  write.csv(core, "output/tables/TableS17b_core_signature.csv", row.names = FALSE)
  log_msg("跨队列核心稳健签名(n≥3队列同向 & FDR<0.05): ", nrow(core), " 个基因")

  save_fig("Fig12_metaZ_volcano", 7, 5.8, {
    print(ggplot(meta, aes(z_meta, -log10(p_meta),
                           color = dir_consistent & fdr < 0.05)) +
          geom_point(size = .6, alpha = .5) +
          scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey75")) +
          labs(x = "Stouffer meta-Z", y = expression(-log[10]~p[meta]),
               color = "Core signature",
               title = paste0("Cross-cohort meta-analysis (", length(zl), " cohorts)"),
               subtitle = paste0("core convergent signature: n = ", nrow(core),
                                 " genes (TableS17b)")) +
          theme(legend.position = "top"))
  })
  if (nrow(core) >= 10) {
    show_genes <- core %>% slice_head(n = 40) %>% pull(gene)
    sign_mat <- zs %>% filter(gene %in% show_genes) %>%
      mutate(s = sign(z)) %>%
      select(gene, cohort, s) %>%
      pivot_wider(names_from = cohort, values_from = s, values_fill = 0) %>%
      column_to_rownames("gene") %>% as.matrix()
    sign_mat <- sign_mat[show_genes, names(zl), drop = FALSE]
    save_fig("FigS24_metaZ_cohort_heatmap", 5.5, 7, {
      print(Heatmap(sign_mat, name = "direction", col = c("-1" = "#00A087",
                    "0" = "grey85", "1" = "#E64B35"),
                    cluster_columns = FALSE, row_names_gp = gpar(fontsize = 7),
                    column_title = "Direction by cohort (top core genes)"))
    })
  }

  quad <- function(lfc1, lfc2, topn = 200) {
    nm1 <- names(lfc1); nm2 <- names(lfc2)
    bg <- length(union(nm1, nm2))
    up1 <- nm1[order(-lfc1)][1:topn]; dn1 <- nm1[order(lfc1)][1:topn]
    up2 <- nm2[order(-lfc2)][1:topn]; dn2 <- nm2[order(lfc2)][1:topn]
    a <- length(intersect(up1, up2)); b <- length(intersect(dn1, dn2))
    ph <- function(k) phyper(k - 1, topn, bg - topn, topn, lower.tail = FALSE)
    tibble(up_up = a, down_down = b, bg = bg,
           hypergeo_p_up = ph(a), hypergeo_p_down = ph(b))
  }
  cv <- safe_csv("output/tables/Table4_callus_validation.csv")
  if (!is.null(cv) && nrow(cv) > 200) {
    q <- quad(setNames(cv$lfc_peri, cv$gene),
              setNames(cv$lfc_callus, cv$gene), topn = 200)
    log_msg("骨痂top200象限富集: up-up = ", q$up_up, " (p=", signif(q$hypergeo_p_up,3),
            "), down-down = ", q$down_down, " (p=", signif(q$hypergeo_p_down,3), ")")
    write.csv(q, "output/tables/TableS17c_callus_quadrant.csv", row.names = FALSE)
    cv2 <- cv %>% mutate(quadrant = case_when(
      lfc_peri > 0 & lfc_callus > 0 ~ "Up-Up", lfc_peri < 0 & lfc_callus < 0 ~ "Down-Down",
      TRUE ~ "Discordant"))
    save_fig("FigS25_callus_quadrant", 6.5, 5, {
      print(ggplot(cv2, aes(lfc_peri, lfc_callus, color = quadrant)) +
            geom_point(size = 1, alpha = .45) +
            scale_color_manual(values = c("Up-Up" = "#E64B35",
                "Down-Down" = "#00A087", "Discordant" = "grey80")) +
            geom_hline(yintercept = 0, color = "grey70") +
            geom_vline(xintercept = 0, color = "grey70") +
            labs(title = "Tail-focused replication (top-200 quadrants)",
                 subtitle = paste0("Up-Up n=", q$up_up, " p=", signif(q$hypergeo_p_up,2),
                                   " | Down-Down n=", q$down_down, " p=",
                                   signif(q$hypergeo_p_down,2)),
                 x = "Periosteum aging log2FC", y = "Callus aging log2FC", color = NULL))
    })
  }
})

## ============================================================
## §7 GSE198666 骨痂时间序列（若有≥2个时间点；v1.1: 兼容 "day 7" 数字在后）
## ============================================================
run_step("timeseries_198", {
  f1 <- "output/tables/TableS0_meta198.csv"; f2 <- "output/seu198.rds"
  if (!all(file.exists(c(f1, f2)))) {
    log_msg("缺 meta198/seu198（主脚本 v3.8 补丁4）→ 跳过时间序列"); return(NULL) }
  meta198 <- read.csv(f1, stringsAsFactors = FALSE)
  txt <- tolower(paste(meta198$title, meta198$source))
  tp <- suppressWarnings(as.numeric(
    str_match(txt, "(\\d+)\\s*(?:d|day)(?:[^a-z]|$)")[, 2]))
  tp <- ifelse(is.na(tp), suppressWarnings(as.numeric(
    str_match(txt, "(?:^|\\s)(?:d|day)\\s*(\\d+)")[, 2])), tp)
  if (n_distinct(na.omit(tp)) < 2) {
    log_msg("GSE198666 仅单一时间点（解析到: ",
            paste(unique(na.omit(tp)), collapse = ","),
            "）→ 跳过衰减曲线（若为单时间点设计属正常）"); return(NULL) }
  meta198$timepoint <- tp
  seu198 <- readRDS(f2)
  seu198$timepoint <- meta198$timepoint[match(seu198$donor_id, meta198$gsm)]
  seu198 <- seu198[, !is.na(seu198$timepoint)]
  prog <- safe_csv("output/tables/Table3_high_confidence_genes.csv")
  prog_genes <- if (!is.null(prog)) intersect(prog$gene, rownames(seu198)) else
    intersect(c("Postn","Col1a1","Sp7","Ibsp","Bglap","Spp1","Mki67","Pdgfra","Dpt"),
              rownames(seu198))
  if (length(prog_genes) < 5) { log_msg("响应程序基因不足"); return(NULL) }
  seu198 <- AddModuleScore(seu198, features = list(prog_genes), name = "RespProg")
  ts_df <- seu198@meta.data %>%
    group_by(age_group, timepoint) %>%
    summarise(score = mean(RespProg1, na.rm = TRUE), n = n(), .groups = "drop")
  peak <- ts_df %>% group_by(age_group) %>%
    summarise(peak = max(score), t_peak = timepoint[which.max(score)], .groups = "drop")
  peak_ratio <- peak$peak[peak$age_group == "Young"] /
                peak$peak[peak$age_group == "Aged"]
  write.csv(ts_df, "output/tables/TableS18_timeseries_scores.csv", row.names = FALSE)
  log_msg("响应曲线峰值 Young/Aged = ", round(peak_ratio, 2))
  save_fig("Fig13_response_timeseries", 6.5, 4.5, {
    print(ggplot(ts_df, aes(timepoint, score, color = age_group, group = age_group)) +
          geom_line(linewidth = 1) + geom_point(size = 2.5) +
          scale_color_manual(values = c("Young" = "#00A087", "Aged" = "#E64B35")) +
          labs(x = "Days post-fracture", y = "Young-response program score",
               color = NULL, title = "Abortive activation: response amplitude over time",
               subtitle = paste0("peak fold Young/Aged = ", round(peak_ratio, 2))) +
          theme(legend.position = "top"))
  })
  rm(seu198); gc(reset = TRUE)
})

## ============================================================
## §8 骨膜分子衰老时钟（默认关闭：训练集 n=3 不成立，不用于论文）
## ============================================================
run_step("aging_clock", {
  if (!DO_CLOCK) { log_msg("DO_CLOCK=FALSE（训练集 n=3，时钟结果不应用于论文）→ 跳过"); return(NULL) }
  if (!file.exists("output/pb297.RData")) {
    log_msg("缺 output/pb297.RData → 跳过时钟"); return(NULL) }
  if (!file.exists(PERI_RDS)) { log_msg("缺 peri RDS → 跳过时钟"); return(NULL) }
  if (!requireNamespace("glmnet", quietly = TRUE)) { log_msg("缺 glmnet"); return(NULL) }
  load("output/pb297.RData")
  cts_tr <- GetAssayData(pb297, assay = "RNA", layer = "counts")
  age_tr <- as.numeric(pb297$age_mo)
  log_msg("时钟训练集: ", ncol(cts_tr), " donors, 月龄范围 ",
          range(age_tr, na.rm = TRUE), "（注意: n<10, 仅供内部探索）")

  peri <- readRDS(PERI_RDS)
  sub <- subset(peri, subset = fraction == "CD45neg" & celltype %in% STROMAL)
  donors <- unique(sub$donor_id)
  cts_q <- do.call(cbind, lapply(donors, function(d)
    Matrix::rowSums(GetAssayData(sub, assay = "RNA", layer = "counts")[
      , sub$donor_id == d, drop = FALSE])))
  colnames(cts_q) <- donors
  age_q <- AGE_MONTHS[sub$age_group[match(donors, sub$donor_id)]]
  rm(peri, sub); gc(reset = TRUE)

  common <- intersect(rownames(cts_tr), rownames(cts_q))
  lcpm <- function(m) log2(sweep(m, 2, colSums(m), "/") * 1e6 + 1)
  L_tr <- lcpm(cts_tr[common, , drop = FALSE])
  L_q  <- lcpm(cts_q[common, , drop = FALSE])
  mad_v <- apply(L_tr, 1, mad)
  vg <- names(sort(mad_v, decreasing = TRUE))[1:min(2000, length(mad_v))]
  ctr <- rowMeans(L_tr[vg, , drop = FALSE]); scl <- apply(L_tr[vg, , drop = FALSE], 1, sd)
  scl[scl == 0] <- 1
  X_tr <- t(sweep(sweep(L_tr[vg, , drop = FALSE], 1, ctr), 1, scl, "/"))
  vg_q <- vg[vg %in% rownames(L_q)]
  X_q  <- t(sweep(sweep(L_q[vg_q, , drop = FALSE], 1, ctr[vg_q]), 1, scl[vg_q], "/"))

  set.seed(9)
  fit <- glmnet::cv.glmnet(X_tr, age_tr, alpha = 0.5, nfolds = min(3, nrow(X_tr)))
  pred <- drop(predict(fit, X_q, s = "lambda.min"))
  ck <- cor.test(pred, age_q, method = "pearson")
  log_msg("时钟外推: predicted vs actual r = ", round(unname(ck$estimate), 2),
          " (p = ", signif(ck$p.value, 3), ") — 仅供内部探索, 不用于论文")
  clock_df <- tibble(donor = donors, age_actual = unname(age_q), age_predicted = pred,
                     group = names(age_q))
  write.csv(clock_df, "output/tables/TableS19_clock.csv", row.names = FALSE)
  save_fig("Fig14_aging_clock", 5.5, 4.8, {
    print(ggplot(clock_df, aes(age_actual, age_predicted, color = group)) +
          geom_point(size = 4) +
          geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
          scale_color_manual(values = c("Young" = "#00A087", "Aged" = "#E64B35"),
                             breaks = c("Young", "Aged")) +   # v1.2(G): 图例 Young 左 Aged 右
          labs(title = "Periosteal molecular aging clock (EXPLORATORY, n=3 train)",
               subtitle = paste0("r = ", round(unname(ck$estimate), 2)),
               x = "Chronological age (months)", y = "Predicted age (months)", color = NULL) +
          theme(legend.position = "top"))
  })
  rm(cts_tr, cts_q, L_tr, L_q); gc(reset = TRUE)
})

## ============================================================
## §9 多轨迹方法交叉验证 + tradeSeq 停滞基因（加固 rho=-0.7 与停滞）
## ============================================================
run_step("trajectory_crossval", {
  if (!file.exists("output/monocle3_cds_periosteum.RData")) {
    log_msg("缺 monocle3 cds → 跳过轨迹交叉验证"); return(NULL) }
  load("output/monocle3_cds_periosteum.RData")
  emb <- as.matrix(reducedDims(cds)[["UMAP"]])
  cl  <- as.character(colData(cds)$celltype)
  grp <- as.character(colData(cds)$group)
  pt_mono <- as.numeric(colData(cds)$pseudotime)

  res_list <- list(Monocle3 = pt_mono)
  if (requireNamespace("slingshot", quietly = TRUE)) {
    sl <- tryCatch(slingshot::slingshot(emb, clusterLabels = cl,
                                        start.clus = "pSSPC_q", approx_points = 200),
                   error = function(e) { log_msg("slingshot失败: ", e$message); NULL })
    if (!is.null(sl)) res_list$Slingshot <- slingshot::slingPseudotime(sl)[, 1]
  } else log_msg("缺 slingshot 包 → 仅 Monocle3")
  if (requireNamespace("destiny", quietly = TRUE)) {
    pt_dpt <- tryCatch({
      dm <- destiny::DiffusionMap(emb, verbose = FALSE)
      tips <- which(cl == "pSSPC_q")[1:min(3, sum(cl == "pSSPC_q"))]
      dpt <- destiny::DPT(dm, tips = tips)
      as.numeric(tryCatch(dpt$dpt, error = function(e) as.matrix(dpt)[, 1]))
    }, error = function(e) { log_msg("DPT失败: ", e$message); NULL })
    if (!is.null(pt_dpt)) res_list$DPT <- pt_dpt
  } else log_msg("缺 destiny 包 → 仅 Monocle3")

  summ <- map_dfr(names(res_list), function(mnm) {
    ptv <- res_list[[mnm]]
    tibble(method = mnm, pt = as.numeric(ptv),
           group = grp[seq_along(ptv)]) %>%
      filter(is.finite(pt)) %>%
      group_by(method, group) %>% summarise(med = median(pt), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = med) %>%
      mutate(delta_young = Young_Fracture - Young_Intact,
             delta_aged  = Aged_Fracture  - Aged_Intact,
             attenuation = delta_young - delta_aged)
  })
  print(as.data.frame(summ))
  write.csv(summ, "output/tables/TableS20_trajectory_methods.csv", row.names = FALSE)
  log_msg("轨迹交叉验证: ", nrow(summ), " 种方法; 老龄推进减弱方向一致 = ",
          all(summ$attenuation < 0) || all(summ$attenuation > 0))

  wilc <- map_dfr(names(res_list), function(mnm) {
    ptv <- res_list[[mnm]]; df <- tibble(pt = as.numeric(ptv), group = grp)
    tt <- wilcox.test(pt ~ group, data = df %>%
      filter(group %in% c("Young_Fracture", "Aged_Fracture")))
    tibble(method = mnm, p_fracture_y_vs_a = signif(tt$p.value, 3))
  })
  write.csv(wilc, "output/tables/TableS20b_trajectory_wilcox.csv", row.names = FALSE)

  save_fig("Fig15_trajectory_crossval", 11, 4.2 * length(res_list), {
    plist <- map(names(res_list), function(mnm) {
      df <- tibble(pt = as.numeric(res_list[[mnm]]), group = grp) %>%
        mutate(group = factor(group, levels = c("Young_Intact","Young_Fracture",
                                                "Aged_Intact","Aged_Fracture")))
      ggplot(df, aes(pt, fill = group)) + geom_density(alpha = .55) +
        scale_fill_manual(values = c("Young_Intact" = "#F39B7F",
          "Young_Fracture" = "#7CAE00", "Aged_Intact" = "#00BFC4",
          "Aged_Fracture" = "#C77CFF")) +
        labs(title = mnm, x = "Pseudotime", y = "Density", fill = NULL)
    })
    print(wrap_plots(plist, ncol = 1) +
          plot_annotation(title = "Aged attenuation of regenerative progression across methods",
                          theme = theme(plot.title = element_text(face = "bold"))))
  })

  if (DO_TRADESEQ && requireNamespace("tradeSeq", quietly = TRUE) &&
      "Slingshot" %in% names(res_list)) {
    counts_tr <- assay(cds, "counts")
    vg2 <- names(sort(apply(log1p(counts_tr), 1, var), decreasing = TRUE))[1:3000]
    cw <- matrix(1, ncol(counts_tr), 1)
    gam <- tradeSeq::fitGAM(counts = counts_tr[vg2, ], pseudotime = res_list$Slingshot,
                            cellWeights = cw, nknots = 6, verbose = FALSE,
                            conditions = factor(grp))
    at <- tradeSeq::conditionTest(gam, lineages = 1)
    stalled <- as.data.frame(at) %>% rownames_to_column("gene") %>%
      filter(pvalue < 0.01) %>% arrange(pvalue)
    write.csv(stalled, "output/tables/TableS21_stalled_genes_tradeSeq.csv",
              row.names = FALSE)
    log_msg("tradeSeq 条件依赖轨迹基因: ", nrow(stalled))
    egost <- enrichGO(stalled$gene[1:min(200, nrow(stalled))], OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL", ont = "BP")
    write.csv(egost@result, "output/tables/TableS21b_stalled_GO.csv", row.names = FALSE)
  }
  rm(cds); gc(reset = TRUE)
})

## ============================================================
## §10 CellChat 差异通讯统计 + 配体活性分析（v1.1: 手动差值热图修复）
## ============================================================
run_step("cellchat_upgrade", {
  if (!file.exists("output/CellChat_periosteum.RData")) {
    log_msg("缺 CellChat RData → 跳过"); return(NULL) }
  load("output/CellChat_periosteum.RData")
  object.list <- list(Young = cc_Y, Aged = cc_A)
  merged.cc <- CellChat::mergeCellChat(object.list, add.names = names(object.list))
  g3 <- tryCatch(CellChat::rankNet(merged.cc, mode = "comparison", stacked = TRUE,
                                   do.stat = TRUE),
                 error = function(e) { log_msg("rankNet失败: ", e$message); NULL })
  ## v1.4(Q): rankNet 空白行（如 CADM/CD226）= 两条龄组信息流均≈0/NA 的并集通路,
  ##          属正常现象; 图中过滤（完整列表见 TableS8/S9 与 CellChat 对象）,
  ##          图高随保留通路数动态调整。
  if (!is.null(g3)) {
    rn_data <- g3$data
    nm_col  <- intersect(c("name", "pathway", "y"), names(rn_data))[1]
    val_col <- intersect(c("contribution", "value", "relative", "count"),
                         names(rn_data))[1]
    if (!is.null(nm_col) && !is.null(val_col)) {
      keep <- rn_data %>%
        filter(!is.na(.data[[val_col]]), .data[[val_col]] > 0) %>%
        pull(!!nm_col) %>% unique()
      drop_lv <- setdiff(unique(rn_data[[nm_col]]), keep)
      if (length(drop_lv) > 0)
        log_msg("rankNet 过滤零/NA 信息流通路（空白行, 正常）: ",
                paste(drop_lv, collapse = ", "))
      rn_data <- rn_data %>% filter(.data[[nm_col]] %in% keep)
      orig_lv <- unique(g3$data[[nm_col]])
      rn_data[[nm_col]] <- factor(rn_data[[nm_col]],
                                  levels = orig_lv[orig_lv %in% keep])
      g3$data <- rn_data
    }
    n_path <- length(unique(g3$data[[nm_col]]))
    g3b <- g3 + theme(axis.text.y = element_text(size = 6),
                      axis.text.x = element_text(size = 8),
                      plot.margin = margin(6, 12, 6, 6))
    save_fig("FigS26_ranknet", 8.5, max(8, 0.2 * n_path + 2), print(g3b))
  }

  ## v1.1: 手动 count/weight 差值热图（规避 netVisual_diffInteraction 的
  ##       "replacement has length zero"——某些细胞类型在单一年龄缺失所致）
  diff_df <- function(mY, mA) {
    ct <- union(rownames(mY), rownames(mA))
    pad <- function(m) { m <- m[ct, ct, drop = FALSE]; m[is.na(m)] <- 0; m }
    as.data.frame(pad(mA) - pad(mY)) %>% rownames_to_column("source") %>%
      pivot_longer(-source, names_to = "target", values_to = "diff")
  }
  d_count  <- diff_df(cc_Y@net$count, cc_A@net$count)
  d_weight <- diff_df(cc_Y@net$weight, cc_A@net$weight)
  p_dc <- ggplot(d_count, aes(target, source, fill = diff)) +
    geom_tile() + scale_fill_gradient2(low = "#00A087", high = "#E64B35") +
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 6)) +
    labs(title = "Interaction count (Aged − Young)", x = NULL, y = NULL, fill = "Δ")
  p_dw <- ggplot(d_weight, aes(target, source, fill = diff)) +
    geom_tile() + scale_fill_gradient2(low = "#00A087", high = "#E64B35") +
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 6)) +
    labs(title = "Interaction weight (Aged − Young)", x = NULL, y = NULL, fill = "Δ")
  save_fig("Fig16_cellchat_diff", 12, 5.5, print(wrap_plots(p_dc, p_dw, ncol = 2)))
  rm(cc_Y, cc_A, merged.cc); gc(reset = TRUE)

  ## v1.1: 配体活性——显式拆分 donor/celltype + 均值比 + 空值守卫
  if (!file.exists(PERI_RDS)) return(NULL)
  peri <- readRDS(PERI_RDS)
  data("CellChatDB.mouse", package = "CellChat")
  lig_all <- unique(na.omit(CellChatDB.mouse$interaction$ligand))
  rec_all <- unique(unlist(strsplit(unique(na.omit(CellChatDB.mouse$interaction$receptor)),
                                    "_")))
  gene_use <- intersect(union(lig_all, rec_all), rownames(peri))
  set.seed(11)
  cs <- colnames(peri)[sample(ncol(peri), min(30000, ncol(peri)))]
  obj <- subset(peri, cells = cs)
  avg <- AverageExpression(obj, features = gene_use,
                           group.by = c("donor_id", "celltype"))$RNA
  long <- avg %>% as.data.frame() %>% rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "dc", values_to = "expr") %>%
    mutate(donor_id = sub("_.*$", "", dc),
           celltype = sub("^[^_]*_", "", dc)) %>%
    left_join(distinct(obj@meta.data, donor_id, age_group), by = "donor_id") %>%
    filter(!is.na(age_group))
  lig_fold <- long %>%
    filter(gene %in% intersect(lig_all, gene_use), celltype %in% SENDERS) %>%
    group_by(celltype, gene) %>%
    summarise(aged_expr  = mean(expr[age_group == "Aged"],  na.rm = TRUE),
              young_expr = mean(expr[age_group == "Young"], na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(aged_expr), is.finite(young_expr), young_expr > 1e-6) %>%
    mutate(fold = aged_expr / young_expr) %>%
    group_by(gene) %>%
    summarise(max_sender_fold = max(fold),
              best_sender = celltype[which.max(fold)], .groups = "drop")
  psub <- subset(obj, subset = celltype %in% c("pSSPC_fib","pSSPC_osteo","pSSPC_q") &
                          age_group == "Aged")
  cnt <- GetAssayData(psub, assay = "RNA", layer = "counts")
  rr <- intersect(rec_all, rownames(cnt))
  rec_stat <- tibble(gene = rr,
                     rec_pct  = Matrix::rowMeans(cnt[rr, , drop = FALSE] > 0) * 100,
                     rec_mean = Matrix::rowMeans(cnt[rr, , drop = FALSE]))
  lr_score <- CellChatDB.mouse$interaction %>%
    select(ligand, receptor) %>% distinct() %>%
    mutate(receptor = as.character(receptor)) %>%
    inner_join(lig_fold, by = c("ligand" = "gene")) %>%
    inner_join(rec_stat, by = c("receptor" = "gene")) %>%
    mutate(score = max_sender_fold * sqrt(rec_pct)) %>% arrange(desc(score))
  write.csv(lr_score, "output/tables/TableS22_ligand_activity.csv", row.names = FALSE)
  log_msg("配体活性 top: "); print(as.data.frame(head(lr_score, 12)))
  ## v1.3(N): 图内仅展示 rec_pct>=5% 的 LR 对（低比例受体噪声大，保留在 TableS22 全表）
  lr_plot <- lr_score %>% filter(rec_pct >= 5)
  save_fig("Fig17_ligand_activity", 8, 6, {
    print(ggplot(head(lr_plot, 20),
                 aes(score, reorder(paste(ligand, receptor, sep = "→"), score),
                     fill = best_sender)) +
          geom_col() + labs(x = "Sender fold × sqrt(receptor pct in aged pSSPC)",
                            y = NULL, fill = "Best sender",
                            title = "Candidate drivers of the aged pSSPC state",
                            subtitle = "simplified NicheNet-style scoring (exploratory)"))
  })
  rm(peri, obj); gc(reset = TRUE)
})

## ============================================================
## §11 coloc.susie 多信号共定位 + 敏感性分析（补救 coloc.abf 全阴）
## ============================================================
run_step("coloc_susie", {
  if (!requireNamespace("coloc", quietly = TRUE)) { log_msg("缺 coloc"); return(NULL) }
  co <- safe_csv("output/tables/Table7_MR_causal_genes.csv")
  if (is.null(co) || nrow(co) == 0) { log_msg("无 MR 阳性基因 → 跳过"); return(NULL) }
  coloc_genes <- unique(co$gene)
  if (length(coloc_genes) == 0) { log_msg("跳过"); return(NULL) }
  if (!requireNamespace("ieugwasr", quietly = TRUE)) { log_msg("缺 ieugwasr"); return(NULL) }
  OUTCOME_ID <- "ebi-a-GCST90038703"; N_GWAS <- 484598; S_GWAS <- 8844/484598
  EqtlN <- c(Fibroblasts = 483, WholeBlood = 670)

  run_one <- function(gene, tissue_name) {
    f <- sprintf("coloc_data/eQTL_%s_%s.csv", gene, tissue_name)
    if (!file.exists(f)) return(NULL)
    eq <- read.csv(f, stringsAsFactors = FALSE) %>%
      mutate(se_e = abs(beta) / pmax(abs(qnorm(pval / 2)), 0.5)) %>%
      filter(is.finite(se_e), se_e > 0)
    gw <- tryCatch(ieugwasr::associations(variants = eq$snp, id = OUTCOME_ID, proxies = 0),
                   error = function(e) NULL)
    if (is.null(gw) || nrow(gw) < 50) return(NULL)
    rsid_col <- intersect(c("rsid","variant","snp","name"), names(gw))[1]
    ea_col <- intersect(c("effect_allele","ea"), names(gw))[1]
    oa_col <- intersect(c("other_allele","nea","oa"), names(gw))[1]
    locus <- gw %>%
      transmute(snp = .data[[rsid_col]], beta_g = beta, se_g = se,
                maf_g = eaf, gw_eff = .data[[ea_col]], gw_oth = .data[[oa_col]]) %>%
      inner_join(eq %>% select(snp, beta_e = beta, se_e, eq_ref = ref, eq_alt = alt),
                 by = "snp") %>%
      mutate(flip = case_when(
        toupper(eq_alt) == toupper(gw_eff) & toupper(eq_ref) == toupper(gw_oth) ~ 1L,
        toupper(eq_alt) == toupper(gw_oth) & toupper(eq_ref) == toupper(gw_eff) ~ -1L,
        TRUE ~ 0L)) %>%
      filter(flip != 0L) %>% mutate(beta_e = beta_e * flip, MAF = coalesce(maf_g, 0.3))
    if (nrow(locus) < 50) return(NULL)
    ## v1.6(T): runsusie 需要 in-sample LD；尝试从 OpenGWAS 参考面板获取，失败则 susie 记 NA
    LD_mat <- tryCatch(ieugwasr::ld_matrix(locus$snp, pop = "EUR"),
                       error = function(e) {
                         log_msg("  LD获取失败(", gene, ",", tissue_name, "): ", e$message); NULL })
    if (!is.null(LD_mat)) {
      cs <- intersect(locus$snp, rownames(LD_mat))
      if (length(cs) >= 50) {
        LD_mat <- LD_mat[cs, cs, drop = FALSE]
        locus <- locus[match(cs, locus$snp), , drop = FALSE]
        log_msg("  LD 矩阵可用: ", length(cs), " SNPs → 尝试 susie")
      } else LD_mat <- NULL
    }
    D1 <- list(beta = locus$beta_e, varbeta = locus$se_e^2, N = EqtlN[[tissue_name]],
               type = "quant", snp = locus$snp, MAF = locus$MAF)
    D2 <- list(beta = locus$beta_g, varbeta = locus$se_g^2, N = N_GWAS, type = "cc",
               snp = locus$snp, MAF = locus$MAF, s = S_GWAS)
    if (!is.null(LD_mat)) { D1$LD <- LD_mat; D2$LD <- LD_mat }
    abf <- tryCatch(coloc::coloc.abf(D1, D2)$summary, error = function(e) NULL)
    ## v1.3(O): 失败原因写入日志（解释 susie 为何全 NA）
    S1 <- tryCatch(coloc::runsusie(D1), error = function(e) {
      log_msg("  runsusie eQTL(", gene, ",", tissue_name, "): ", e$message); NULL })
    S2 <- tryCatch(coloc::runsusie(D2), error = function(e) {
      log_msg("  runsusie GWAS(", gene, ",", tissue_name, "): ", e$message); NULL })
    sus_res <- NULL
    if (!is.null(S1) && !is.null(S2) && length(S1@pip) > 0 && length(S2@pip) > 0) {
      sus_res <- tryCatch(coloc::coloc.susie(S1, S2)$summary, error = function(e) NULL)
    }
    tibble(gene = gene, tissue = tissue_name, nsnps = nrow(locus),
           PP_H4_abf = if (!is.null(abf)) unname(abf["PP.H4.abf"]) else NA_real_,
           PP_H4_susie = if (!is.null(sus_res)) unname(sus_res["PP.H4.abf"]) else NA_real_,
           n_signals_eqtl = if (!is.null(S1)) length(S1@pip[S1@pip > 0.5]) else NA_integer_,
           n_signals_gwas = if (!is.null(S2)) length(S2@pip[S2@pip > 0.5]) else NA_integer_,
           susie_pass = !is.null(sus_res) && unname(sus_res["PP.H4.abf"]) > 0.75)
  }
  cs_res <- map_dfr(coloc_genes, function(g)
    map_dfr(c("Fibroblasts", "WholeBlood"), function(t)
      tryCatch(run_one(g, t), error = function(e) {
        log_msg("coloc.susie 失败(", g, ",", t, "): ", e$message); NULL })))
  write.csv(cs_res, "output/tables/TableS23_coloc_susie.csv", row.names = FALSE)
  print(as.data.frame(cs_res))
  log_msg("coloc.susie 通过(PP.H4>0.75): ", sum(cs_res$susie_pass, na.rm = TRUE),
          "；abf→susie PP.H4 提升的基因数: ",
          sum(cs_res$PP_H4_susie > cs_res$PP_H4_abf, na.rm = TRUE))
  ## v1.2(J): susie NA 补 "n.e." 标注; 图例保留双色; 防裁剪 + 动态图高
  cs_plot <- cs_res %>%
    pivot_longer(c(PP_H4_abf, PP_H4_susie), names_to = "method", values_to = "PPH4") %>%
    mutate(lab = paste(gene, tissue),
           method = factor(method, levels = c("PP_H4_abf", "PP_H4_susie"),
                           labels = c("coloc.abf", "coloc.susie")))
  sus_na <- cs_res %>% filter(is.na(PP_H4_susie)) %>% mutate(lab = paste(gene, tissue))
  save_fig("FigS27_coloc_susie", 8.5, max(4.5, 0.7 * nrow(cs_res) + 2.5), {
    print(ggplot(cs_plot, aes(PPH4, reorder(lab, PPH4), color = method)) +
          geom_vline(xintercept = 0.75, linetype = 2, color = "grey50") +
          geom_point(size = 3, position = position_dodge(width = 0.5), na.rm = TRUE) +
          { if (nrow(sus_na) > 0)
              geom_text(data = sus_na, aes(x = 0.055, y = lab, label = "susie: n.e."),
                        inherit.aes = FALSE, hjust = 0, size = 2.8, color = "grey45") } +
          scale_color_manual(values = c("coloc.abf" = "#4DBBD5", "coloc.susie" = "#E64B35"),
                             drop = FALSE) +
          scale_y_discrete(expand = expansion(mult = c(0.12, 0.12))) +
          coord_cartesian(clip = "off") +
          labs(x = "PP.H4", y = NULL, color = NULL,
               title = "coloc.abf vs coloc.susie (multi-signal)",
               caption = "n.e. = not estimable (in-sample LD unavailable → susie not run; single-variant coloc.abf reported)") +
          theme(legend.position = "top",
                plot.margin = margin(12, 16, 6, 6)))
  })
})

## ============================================================
## §12 汇总：写给 cover letter 的定量资产清单
## ============================================================
run_step("summary_assets", {
  grab <- function(path, note) if (file.exists(path))
    tibble(file = path, note = note) else NULL
  assets <- bind_rows(
    grab("output/tables/TableS12_stability_bootstrap.csv", "方向敏感性（非独立重复验证）"),
    grab("output/tables/TableS13_empirical_pvalues.csv", "置换/精确检验 p 值"),
    grab("output/tables/TableS15_multisig_donor.csv", "多衰老签名 concordance"),
    grab("output/tables/TableS16_pathway_GSEA.csv", "通路级 GSEA"),
    grab("output/tables/TableS17_metaZ.csv", "跨队列 meta-Z"),
    grab("output/tables/TableS17b_core_signature.csv", "核心稳健签名"),
    grab("output/tables/TableS17c_callus_quadrant.csv", "骨痂尾部象限富集"),
    grab("output/tables/TableS18_timeseries_scores.csv", "响应衰减曲线"),
    grab("output/tables/TableS20_trajectory_methods.csv", "多轨迹交叉验证"),
    grab("output/tables/TableS22_ligand_activity.csv", "配体活性"),
    grab("output/tables/TableS23_coloc_susie.csv", "coloc.susie"))
  if (!is.null(assets)) {
    write.csv(assets, "output/tables/TableS24_new_assets_index.csv", row.names = FALSE)
    log_msg("新增定量资产 ", nrow(assets), " 项（TableS24 索引）")
  }
  log_msg("========== strengthening v1.1 流程结束 ==========")
})
