#!/usr/bin/env Rscript
## ============================================================
## run_periosteum_aging19.R v3.8 —— 单细胞与bulk转录组揭示骨膜衰老的
##   趋同衰老重塑与再生轨迹停滞（图表修正版 + 强化脚本配套修正版）
## 发现层：GSE280914（King et al., Bone 2025）
## 设计：年轻(3-4mo) vs 老龄(20-24mo) × 完整 vs 骨折(d3) × CD45+/CD45- scRNA
##
## v3.8 变更（相对 v3.7c-figures）：
##  A) §8b GSE198666 末尾新增补丁4存档（TableS0_meta198.csv + seu198.rds）；
##  B) §9 GSE297256 重写：主路径按解包目录名硬解析 Y/M/O + 直接 ReadMtx + 免 harmony
##     （修复 v3.7 实测 "invalid 'data'"），保留原 title 解析为回退路径；
##     末尾新增补丁3存档 pb297.RData；
##  C) §9 GSE232516 / §9b GSE278165：babelgene::orthologs 失败时回退 toupper 映射
##     （修复 "no orthologs found or the genes are not valid human genes"）；
##  D) §12 MR FigS10：position_dodge2(height=0.62) → width=0.62（height 非该函数参数，
##     是 v3.7 实测报错根源），并加 position_dodge 兜底；
##  E) 第14节末尾新增补丁1/2存档（peri_seurat.rds + TableS0_full_DE_stats.csv），
##     为 run_periosteum_strengthening.R 提供输入；均带 exists() 防呆。
## 其余内容同 v3.7c（v3.7c 修复清单见下方注释）。
##
## ---- v3.7c 修复清单（继承） ----
##  1) Fig1b 空白：显式 group + position_dodge2(preserve="single")，并 stopifnot 非空；
##  2) GSEA 统一 pvalueCutoff=1 + eps=1e-10；新增 gsea_ok() 守卫出图；
##  3) donor n=1v1 不再展示 AUC=1；供体面板按 n 自动选择 boxplot/jitter 或 point；
##  4) 全部含年龄分组的图统一 Young→Aged 顺序；
##  5) Fig3b/Fig3c 标题/副标题截断修复；
##  6) Fig9/FigS10：recode 硬映射结局短标签 + 映射校验 + 图例两行置底；
##  7) Fig6b 弦图回退参数并把真实报错写入日志；
##  8) FigS4 facet 按 Young→Aged 排序；
##  9) 拟时 Fig7：root 星标 + 数字里程碑与箭头；
## 10) GSE297256 文件名 Y/M/O 标签回退；GSE232516 design 变量还原；
## 11) GSE278165 人参考自动聚类 + TransferData 显式 weight.reduction="cca"。
## ============================================================

## ---------------- 0. 环境 ----------------
suppressPackageStartupMessages({
  library(GEOquery); library(Seurat); library(DESeq2); library(limma)
  library(Matrix); library(tidyverse); library(data.table)
  library(patchwork); library(ComplexHeatmap); library(EnhancedVolcano)
  library(clusterProfiler); library(org.Mm.eg.db); library(org.Hs.eg.db)
  library(Biobase); library(cowplot)
})
library(enrichplot); library(pROC)
if (!requireNamespace("harmony", quietly = TRUE)) install.packages("harmony", type = "binary")
library(harmony)
if (!requireNamespace("babelgene", quietly = TRUE)) install.packages("babelgene")
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; count <- dplyr::count; summarise <- dplyr::summarise
arrange <- dplyr::arrange; desc <- dplyr::desc

for (d in c("data/GSE280914","data/GSE297256","data/GSE198666","data/GSE232516",
            "data/GSE278165","output/figures","output/tables","output/logs","coloc_data"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
log_file <- sprintf("output/logs/run_periosteum_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste(..., collapse = " "))
  cat(msg, "\n"); write(msg, log_file, append = TRUE)
}
run_step <- function(nm, expr) {
  log_msg("---- STEP: ", nm)
  tryCatch(eval(substitute(expr), envir = parent.frame()),
           error = function(e) { log_msg("!! ERROR in ", nm, ": ", e$message); NULL })
}
save_fig <- function(name, width, height, expr, dpi = 600, dir = "output/figures") {
  e <- substitute(expr); envir <- parent.frame()
  pdf(file.path(dir, paste0(name, ".pdf")), width = width, height = height)
  print(eval(e, envir)); dev.off()
  tiff(file.path(dir, paste0(name, ".tiff")), width = width, height = height,
       units = "in", res = dpi, compression = "lzw", type = "cairo")
  print(eval(e, envir)); dev.off()
  log_msg("图已生成(PDF+600dpi TIFF): ", name)
}
dl_with_retry <- function(url, dest, tries = 6, timeout_sec = 1800) {
  if (file.exists(dest) && file.size(dest) > 1000) {
    log_msg("已存在，跳过: ", basename(dest)); return(invisible(TRUE))
  }
  url <- sub("^ftp://", "https://", url)
  have_curl <- requireNamespace("curl", quietly = TRUE)
  if (!have_curl) log_msg("!! 建议 install.packages('curl') 以启用断点续传")
  for (i in seq_len(tries)) {
    ok <- tryCatch({
      if (have_curl) {
        h <- curl::new_handle()
        curl::handle_setopt(h, timeout = timeout_sec, connecttimeout = 60,
                            resume_from = ifelse(file.exists(dest), file.size(dest), 0))
        curl::curl_download(url, dest, mode = "wb", handle = h)
      } else {
        download.file(url, dest, method = "libcurl", mode = "wb", timeout = timeout_sec)
      }
      TRUE
    }, error = function(e) { log_msg("  错误(第", i, "次): ", e$message); FALSE })
    if (ok && file.exists(dest) && file.size(dest) > 1000) {
      log_msg("  完成: ", basename(dest), " (", round(file.size(dest)/1e6, 1), " MB)")
      return(invisible(TRUE))
    }
    Sys.sleep(5)
  }
  log_msg("!! 下载失败（可浏览器手动下载后放到: ", dest, "）"); invisible(FALSE)
}
options(timeout = 1800)

## ---- 通用小工具（v3.7c/v3.8） ----
set_pb_meta <- function(pb, meta, by = "donor_id") {
  md <- data.frame(cell = colnames(pb), stringsAsFactors = FALSE)
  names(md)[1] <- by
  md <- md %>% left_join(meta, by = by)
  rownames(md) <- md[[by]]
  pb@meta.data <- md
  pb
}
gsea_ok <- function(g) {
  !is.null(g) && inherits(g, "gseaResult") && nrow(g@result) > 0 && !is.na(g@result$NES[1])
}
mk_donor_panel <- function(df, ylab = "Mean SenMayo score", title = "") {
  df$age_group <- factor(df$age_group, levels = c("Young","Aged"))
  n_by <- table(df$age_group)
  if (all(n_by >= 2)) {
    ggplot(df, aes(age_group, SenMayo, fill = age_group)) +
      geom_boxplot() + geom_jitter(width = .08, size = 2.5) +
      scale_fill_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
      theme(legend.position = "none") +
      labs(title = title, x = NULL, y = ylab)
  } else {
    ggplot(df, aes(age_group, SenMayo, color = age_group)) +
      geom_point(size = 3) +
      scale_color_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
      theme(legend.position = "none") +
      labs(title = title, x = NULL, y = ylab,
           caption = paste0("donor n = ", paste(n_by, collapse = " vs "), " (descriptive)"))
  }
}
spearman_rho <- function(x, y) {
  r <- suppressWarnings(tryCatch(cor.test(x, y, method = "spearman"),
                                 error = function(e) NULL))
  if (is.null(r)) NA_real_ else unname(r$estimate)
}
## v3.8：babelgene 失败时的安全回退（toupper 直映射，探索性可接受）
safe_orthologs <- function(genes, species = "mouse") {
  og <- tryCatch(babelgene::orthologs(genes = genes, species = species, human = TRUE),
                 error = function(e) {
                   log_msg("babelgene 失败（", e$message, "）→ 回退 toupper 映射"); NULL })
  if (is.null(og) || nrow(og) == 0) return(NULL)
  og
}

log_msg("========== 流程开始 (periosteum_aging v3.8, GSE280914) ==========")

## ---------------- ★全局配置 ----------------
pick_outcomes <- function(keyword = "fracture") {
  library(TwoSampleMR)
  ao <- available_outcomes()
  hit <- ao[str_detect(tolower(ao$trait), keyword), c("id","trait","sample_size","ncase","ncontrol")]
  print(as.data.frame(hit)); invisible(hit)
}
OUTCOME_IDS <- c(
  "ebi-a-GCST90038703",   # 主结局：Fractures（医院诊断） N=484,598, cases=8,844
  "ebi-a-GCST90038705",   # 部位特异：前臂/腕骨折（Colles） cases=2,291
  "ebi-a-GCST006980",     # 敏感性：Fractures（含自报） cases=53,184
  "ebi-a-GCST90029004",   # 定量：跟骨eBMD N=583,314
  "ieu-a-980",            # 定量：股骨颈BMD（GEFOS经典） N=32,735
  "ebi-a-GCST90038656"    # 三级：骨质疏松诊断 cases=7,751
)
OUTCOME_SHORT <- c(
  "ebi-a-GCST90038703" = "Fracture",
  "ebi-a-GCST90038705" = "Forearm fx",
  "ebi-a-GCST006980"   = "Fracture(self-rpt)",
  "ebi-a-GCST90029004" = "eBMD",
  "ieu-a-980"          = "FN-BMD",
  "ebi-a-GCST90038656" = "Osteoporosis"
)
MR_EXTRA_GENES <- c("SPP1","CSF1","TGFB1","IL6","TNF")
DO_ATLAS_VALIDATION <- TRUE
DO_CALLUS_VALIDATION <- TRUE
DO_HUMAN_ATLAS <- TRUE
FRAC_POS <- c("cd45 ?\\+","cd45 ?pos","cd45p","cd45 ?plus","cd45 ?positive")
FRAC_NEG <- c("cd45 ?-","cd45 ?neg","cd45n","cd45 ?minus","cd45 ?negative","cd45 ?deplet")

## ---------------- 1. GEO 通用获取 + 多格式读取 ----------------
fetch_geo <- function(gse, destdir) {
  es <- getGEO(GEO = gse, GSEMatrix = TRUE, getGPL = FALSE)
  eset <- es[[1]]
  pd <- Biobase::pData(eset)
  show_cols <- c("geo_accession","title","source_name_ch1")
  print(as.data.frame(pd[, intersect(show_cols, colnames(pd)), drop = FALSE]))
  ch_cols <- colnames(pd)[str_detect(colnames(pd), "characteristics")]
  if (length(ch_cols))
    print(as.data.frame(pd[, c("geo_accession", ch_cols), drop = FALSE]))
  supp_cols <- colnames(pd)[str_detect(colnames(pd), "supplementary_file")]
  supp <- unique(na.omit(unlist(pd[, supp_cols, drop = FALSE])))
  rawdir <- file.path(destdir, "RAW")
  dir.create(rawdir, showWarnings = FALSE, recursive = TRUE)
  for (u in supp) {
    dest <- file.path(rawdir, basename(u))
    dl_with_retry(u, dest)
    if (str_detect(dest, "\\.tar(\\.gz)?$") && file.exists(dest)) {
      exdir <- paste0(dest, "_untar"); dir.create(exdir, showWarnings = FALSE)
      untar(dest, exdir = exdir)
    }
  }
  log_msg(gse, " 下载完成: ", length(list.files(rawdir, recursive = TRUE)), " 个文件")
  invisible(list(pd = pd, rawdir = rawdir))
}

read_sample_any <- function(rawdir, gsm, m, mt.pattern = "^mt-") {
  f <- list.files(rawdir, pattern = gsm, full.names = TRUE, recursive = TRUE)
  f <- f[!dir.exists(f)]
  if (!length(f)) {  ## tar 解包后文件名常不含 gsm，回退到 <gsm>*.tar.gz_untar 目录
    d <- list.files(rawdir, pattern = paste0(gsm, ".*untar"), full.names = TRUE,
                    recursive = TRUE, include.dirs = TRUE)
    d <- d[dir.exists(d)]
    if (length(d)) f <- list.files(d[1], full.names = TRUE, recursive = TRUE)
  }
  if (!length(f)) stop("找不到样本文件: ", gsm)
  counts <-
    if (any(str_detect(basename(f), "matrix\\.mtx"))) {
      pick <- function(pat) {
        cand <- f[str_detect(basename(f), pat)]
        if (!length(cand)) return(NA_character_)
        filtered <- cand[str_detect(basename(cand), "filtered")]
        if (length(filtered)) return(filtered[1])
        cand[1]
      }
      fb <- pick("barcodes"); ff <- pick("features|genes"); fm <- pick("matrix\\.mtx")
      if (any(is.na(c(fb, ff, fm))))
        stop("样本文件不全: ", gsm, " | 找到: ", paste(basename(f), collapse = ", "))
      ReadMtx(mtx = fm, cells = fb, features = ff, feature.column = 2)
    } else if (any(str_detect(basename(f), "\\.h5$"))) {
      if (!requireNamespace("hdf5r", quietly = TRUE)) install.packages("hdf5r")
      h5f <- f[str_detect(basename(f), "\\.h5$")][1]
      mat <- Read10X_h5(h5f)
      if (is.list(mat)) mat <- mat[[1]]
      mat
    } else if (any(str_detect(basename(f), "\\.rds$"))) {
      return(readRDS(f[str_detect(basename(f), "\\.rds$")][1]))
    } else {
      matf <- f[which.max(file.size(f))]
      log_msg("  读取矩阵文件: ", basename(matf), " (", round(file.size(matf)/1e6, 1), " MB)")
      df <- data.table::fread(matf, data.table = FALSE)
      genes <- df[[1]]
      mat <- as.matrix(df[, -1, drop = FALSE])
      rownames(mat) <- genes
      mat <- mat[!duplicated(genes) & !is.na(genes) & genes != "", , drop = FALSE]
      Matrix(mat, sparse = TRUE)
    }
  obj <- CreateSeuratObject(counts, project = m$sample_id, min.cells = 3, min.features = 200)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt.pattern)
  obj$age_group <- m$age_group
  obj$injury    <- m$injury
  obj$fraction  <- m$fraction
  obj$donor_id  <- m$donor_id
  obj$sample_id <- m$sample_id
  subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
}

## ---------------- 2. GSE280914 元数据 ----------------
fetch_dsc <- run_step("fetch_GSE280914", fetch_geo("GSE280914", "data/GSE280914"))
if (is.null(fetch_dsc)) stop("GSE280914 获取失败")
pd_dsc <- fetch_dsc$pd

pos_pat <- paste(FRAC_POS, collapse = "|")
neg_pat <- paste(FRAC_NEG, collapse = "|")
meta_dsc <- tibble(
  gsm = pd_dsc$geo_accession, title = pd_dsc$title,
  source = pd_dsc$source_name_ch1,
  txt = tolower(paste(title, source))) %>%
  mutate(
    age_group = case_when(
      str_detect(txt, "aged|old|\\b(1[89]|2[0-4]) ?m|geriatric") ~ "Aged",
      str_detect(txt, "young|\\b[2-6] ?m|adult") ~ "Young",
      TRUE ~ NA_character_),
    injury = case_when(
      str_detect(txt, "fractur|fx|post.?injur|day ?3|\\bd3\\b") ~ "Fracture",
      str_detect(txt, "intact|naive|uninjur|sham|contralateral") ~ "Intact",
      TRUE ~ NA_character_),
    fraction = case_when(
      str_detect(txt, "non-?immune") ~ "CD45neg",
      str_detect(txt, "immune") ~ "CD45pos",
      str_detect(txt, pos_pat) ~ "CD45pos",
      str_detect(txt, neg_pat) ~ "CD45neg",
      TRUE ~ "CD45neg"),
    sample_id = gsm, donor_id = gsm)
print(as.data.frame(meta_dsc %>% select(gsm, title, age_group, injury, fraction)))
if (any(is.na(meta_dsc$age_group)) || any(is.na(meta_dsc$injury)))
  stop("元数据解析不完整 → 请override")
stopifnot(all(table(meta_dsc$age_group, meta_dsc$injury) >= 2), nrow(meta_dsc) == 8)
log_msg("GSE280914 分组表: "); print(table(meta_dsc$age_group, meta_dsc$injury, meta_dsc$fraction))

## ---------------- 3. 单细胞预处理与整合 ----------------
objs <- run_step("load_samples", setNames(lapply(seq_len(nrow(meta_dsc)), function(i)
  read_sample_any(fetch_dsc$rawdir, meta_dsc$gsm[i], meta_dsc[i, ], mt.pattern = "^mt-")),
  meta_dsc$sample_id))
if (is.null(objs)) stop("load_samples 失败")
log_msg("样本载入完成，总细胞数: ", sum(sapply(objs, ncol)))

peri <- run_step("integrate", {
  m <- reduce(objs, merge) %>% NormalizeData() %>% FindVariableFeatures() %>%
       ScaleData() %>% RunPCA()
  m <- RunHarmony(m, group.by.vars = "donor_id")
  m %>% FindNeighbors(reduction = "harmony") %>% FindClusters(resolution = 0.5) %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
})
if (is.null(peri)) stop("integrate 失败")
print(table(peri$sample_id)); print(table(Idents(peri)))

## ---------------- 4. 骨膜细胞注释 ----------------
peri_markers <- c(
  "Pdgfra","Col1a1","Col3a1","Prrx1","Ctsk","Postn","Dpt","Itm2a","Lrp1","Anpep",
  "Acta2","Lepr","Gli1","Thy1",
  "Sp7","Runx2","Bglap","Ibsp","Alpl","Sox9","Acan","Col2a1","Col10a1",
  "Ptprc","Cd3e","Cd4","Cd8a","Cd19","Csf1r","Adgre1","S100a8","S100a9","Ly6g","Ncr1",
  "Pecam1","Kdr","Cdh5","Vwf","Rgs5","Cspg4","Pdgfrb",
  "Hbb-b1","Hba-a1","Myog","Myh1","Tubb3","S100b",
  "Mki67","Top2a","Pcna")
peri_markers <- peri_markers[peri_markers %in% rownames(peri)]
save_fig("FigS0_umap", 7, 6, { print(DimPlot(peri, reduction = "umap", label = TRUE)) })
save_fig("FigS1_dotplot_markers", 13.5, 7, {
  print(DotPlot(peri, features = peri_markers) + RotatedAxis() +
          theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8)))
})
log_msg("★ 请查看 FigS1 后确认下方 cluster_names（剔除肌肉/骨髓污染簇）")

cluster_names <- c(
  "0"  = "Neutroph","1" = "Neutroph","2" = "pSSPC_fib","3" = "Neutroph",
  "4"  = "Neutroph","5" = "Macroph","6" = "Neutroph","7" = "Neutroph",
  "8"  = "Cycling","9" = "Macroph","10" = "pSSPC_fib","11" = "Neutroph",
  "12" = "Cycling","13" = "Cycling","14" = "Cycling","15" = "Endo",
  "16" = "Neutroph","17" = "pSSPC_osteo","18" = "Chondro","19" = "Neutroph",
  "20" = "pSSPC_q","21" = "Cycling","22" = "Cycling","23" = "LowQ_removed")
stopifnot(all(names(cluster_names) %in% levels(Idents(peri))))
Idents(peri) <- peri$seurat_clusters
peri <- RenameIdents(peri, cluster_names)
peri$celltype <- as.character(Idents(peri))
peri <- subset(peri, subset = celltype != "LowQ_removed")
peri <- JoinLayers(peri)
log_msg("注释后细胞数: ", ncol(peri))

peri$group4 <- factor(paste(peri$age_group, peri$injury, sep = " "),
                      levels = c("Young Intact","Young Fracture","Aged Intact","Aged Fracture"))
ct_cols <- c("Neutroph"="#E64B35","Macroph"="#F39B7F","pSSPC_q"="#B09C85","pSSPC_fib"="#4DBBD5",
             "pSSPC_osteo"="#3C5488","Chondro"="#91D1C2","Endo"="#00A087","Cycling"="#7E6148")
save_fig("FigS0b_umap_celltype", 8, 6.5, {
print(DimPlot(peri, group.by = "celltype", cols = ct_cols, label = TRUE, repel = TRUE,
              label.size = 4.5, pt.size = 0.25) +
      ggtitle("Periosteal cell atlas — GSE280914") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")))
})
save_fig("FigS0c_umap_by_group", 13.5, 4.2, {
print(DimPlot(peri, group.by = "celltype", cols = ct_cols, split.by = "group4",
              ncol = 4, pt.size = 0.15) +
      ggtitle("Cell type distribution by age × injury") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 10)),
            plot.margin = margin(8, 8, 8, 8),
            strip.text = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 8)))
})
print(table(peri$donor_id, peri$celltype))

## ---------------- 5. 群落级重塑（故事线①） ----------------
comp <- peri@meta.data %>%
  count(donor_id, age_group, injury, fraction, celltype) %>%
  group_by(donor_id, fraction) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
save_fig("Fig1_composition", 9, 5, {
print(ggplot(comp, aes(factor(donor_id), pct, fill = celltype)) +
      geom_col() + coord_flip() +
      facet_wrap(~ fraction + age_group + injury, scales = "free_y") +
      labs(x = NULL, y = "Percentage of cells (%)", fill = NULL))
})
ya <- c("Young","Aged")
fib_d <- comp %>% filter(celltype == "pSSPC_fib", injury == "Intact",
                         fraction == "CD45neg", age_group %in% ya) %>%
  mutate(age_group = factor(age_group, levels = ya))
fib_test <- tryCatch(wilcox.test(pct ~ age_group, data = fib_d), error = function(e) NULL)
imm_pct <- comp %>% filter(fraction == "CD45pos", injury == "Intact", age_group %in% ya) %>%
  group_by(donor_id, age_group) %>% mutate(age_group = factor(age_group, levels = ya)) %>%
  summarise(imm_pct = sum(pct[celltype %in% c("Neutroph","Macroph")]), .groups = "drop")
imm_test <- tryCatch(wilcox.test(imm_pct ~ age_group, data = imm_pct), error = function(e) NULL)
log_msg("CD45neg完整骨膜 pSSPC_fib%（中位 Young/Aged）: ",
        paste(round(tapply(fib_d$pct, fib_d$age_group, median), 1), collapse = " / "),
        if (!is.null(fib_test)) paste0(" | wilcox p = ", signif(fib_test$p.value, 3)) else " | n=1v1无法检验",
        "；CD45pos免疫%: ", paste(round(tapply(imm_pct$imm_pct, imm_pct$age_group, median), 1), collapse = " / "),
        if (!is.null(imm_test)) paste0(" | wilcox p = ", signif(imm_test$p.value, 3)) else " | n=1v1无法检验")
write.csv(comp, "output/tables/TableS1_composition.csv", row.names = FALSE)

fig1b <- comp %>%
  filter(fraction == "CD45neg", injury %in% c("Intact","Fracture"),
         celltype %in% c("pSSPC_q","pSSPC_fib","pSSPC_osteo","Cycling")) %>%
  mutate(celltype = factor(celltype, levels = c("pSSPC_q","pSSPC_fib","pSSPC_osteo","Cycling")),
         age_group = factor(age_group, levels = c("Young","Aged")),
         injury    = factor(injury,    levels = c("Intact","Fracture")))
stopifnot("Fig1b 数据为空，请检查 comp" = nrow(fig1b) > 0)
save_fig("Fig1b_fibrogenic_activation", 9, 4.6, {
print(ggplot(fig1b, aes(celltype, pct, fill = age_group, group = age_group)) +
      geom_col(position = position_dodge2(width = 0.78, preserve = "single"), width = 0.62) +
      geom_text(aes(label = sprintf("%.1f", pct)),
                position = position_dodge2(width = 0.78, preserve = "single"),
                vjust = -0.4, size = 3) +
      facet_wrap(~ injury, nrow = 1) +
      scale_fill_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.22))) +
      labs(x = NULL, y = "Cells (%)", fill = NULL,
           title = "Fibrogenic progenitor activation (CD45neg)",
           caption = "CD45neg fraction; donor n = 1 per age × injury (descriptive)") +
      theme(legend.position = "top",
            plot.title = element_text(size = 11, face = "bold"),
            plot.caption = element_text(size = 7, color = "grey40"),
            strip.text = element_text(face = "bold")))
})

## ---------------- 6. 衰老对比 DE ----------------
STROMAL <- c("pSSPC_q","pSSPC_fib","pSSPC_osteo","Endo","Cycling","Chondro")
DE_SET <- c("pSSPC_q","pSSPC_fib","pSSPC_osteo","Endo","Chondro")
qc_imm <- mean(peri@meta.data$celltype[peri@meta.data$fraction=="CD45neg" &
                                         peri@meta.data$injury=="Intact"] %in% c("Neutroph","Macroph"))
log_msg("QC: CD45neg-intact中免疫簇占比 = ", round(qc_imm*100, 1), "%（v3.3起从间充质分析中剔除）")
peri_intact <- subset(peri, subset = injury == "Intact" & fraction == "CD45neg" &
                        celltype %in% DE_SET)
n_by_age <- peri_intact@meta.data %>% distinct(donor_id, age_group) %>% count(age_group) %>% pull(n)
DE_MODE <- if (all(n_by_age >= 2)) "pseudobulk" else "cell_level"
log_msg("CD45neg完整骨膜每组【文库】数: ", paste(n_by_age, collapse = " / "),
        "（对应 Young/Aged） → DE_MODE = ", DE_MODE)

if (DE_MODE == "pseudobulk") {
  pb <- AggregateExpression(peri_intact, group.by = "donor_id", return.seurat = TRUE)
  mt <- peri_intact@meta.data %>% distinct(donor_id, .keep_all = TRUE)
  pb <- set_pb_meta(pb, mt, by = "donor_id")
  pb$age_group <- factor(pb$age_group, levels = c("Young","Aged"))
  dds <- run_step("DESeq2_aging", {
    DESeqDataSetFromMatrix(GetAssayData(pb, layer = "counts"),
                           pb@meta.data, design = ~ age_group) %>% DESeq()
  })
  if (is.null(dds)) stop("DESeq2_aging 失败")
  res_age <- results(dds, contrast = c("age_group", "Aged", "Young"))
  res_age_df <- as.data.frame(res_age) %>% rownames_to_column("gene")
} else {
  log_msg("!! 每组<2文库 → 细胞级DE（探索性）")
  res_age <- NULL
  Idents(peri_intact) <- peri_intact$age_group
  fm <- FindMarkers(peri_intact, ident.1 = "Aged", ident.2 = "Young",
                    logfc.threshold = 0, min.pct = 0.01, test.use = "wilcox")
  gene_rank <- tryCatch({
    if (!requireNamespace("presto", quietly = TRUE)) stop("presto not installed")
    wauc <- presto::wilcoxauc(peri_intact, group_by = "age_group") %>%
      filter(group == "Aged") %>% select(feature, auc)
    x <- wauc$auc - 0.5; names(x) <- wauc$feature; x
  }, error = function(e) {
    log_msg("  presto不可用 → 退回sign*-log10(p): ", e$message)
    x <- sign(fm$avg_log2FC) * -log10(fm$p_val + 1e-300); names(x) <- rownames(fm); x
  })
  res_age_df <- fm %>% rownames_to_column("gene") %>%
    transmute(gene, log2FoldChange = avg_log2FC, padj = p_val_adj,
              pvalue = p_val, stat = unname(gene_rank[rownames(fm)])) %>%
    filter(!is.na(stat))
}
deg_age <- res_age_df %>% filter(padj < 0.05, abs(log2FoldChange) > 0.25) %>% arrange(padj)
noise_pat <- "^Hist[0-9]|^mt-|^Mt-|^Rpl|^Rps|^Gm[0-9]|^AC[0-9]|^BC[0-9]|^AI[0-9]|^AY[0-9]"
n_before <- nrow(deg_age)
deg_age <- deg_age[!grepl(noise_pat, deg_age$gene), ]
res_age_df <- res_age_df[!grepl(noise_pat, res_age_df$gene), ]
log_msg("噪声基因剔除后 DEG数: ", nrow(deg_age), "（原 ", n_before, "）")
write.csv(deg_age, "output/tables/Table1_aging_DEGs.csv", row.names = FALSE)

volc_df <- res_age_df %>%
  mutate(padj = ifelse(is.na(padj), 1, padj),
         sig = case_when(padj < 0.05 & log2FoldChange >  0.25 ~ "Up in Aged",
                         padj < 0.05 & log2FoldChange < -0.25 ~ "Down in Aged",
                         TRUE ~ "NS"))
save_fig("FigS6_volcano_aging_DE", 7, 6, {
print(ggplot(volc_df, aes(log2FoldChange, -log10(pvalue), color = sig)) +
      geom_point(size = 0.7, alpha = 0.6) +
      scale_color_manual(values = c("Up in Aged"="#E64B35","Down in Aged"="#00A087","NS"="grey80")) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = 2, color = "grey40") +
      labs(x = "log2 fold-change (Aged vs Young)", y = expression(-log[10]~p), color = NULL,
           title = "Periosteal stromal aging DE (CD45neg, intact)") +
      theme(legend.position = "top", plot.title = element_text(size = 11, face = "bold")))
})
hm_genes <- bind_rows(head(deg_age, 25), tail(deg_age, 25)) %>% pull(gene) %>% unique()
hm_avg <- AverageExpression(peri, features = hm_genes, group.by = "sample_id", assays = "RNA")$RNA
hm_mat <- t(scale(t(as.matrix(hm_avg))))
hm_meta <- peri@meta.data %>% distinct(sample_id, age_group, injury, fraction)
ha <- HeatmapAnnotation(
  Age = hm_meta$age_group[match(colnames(hm_mat), hm_meta$sample_id)],
  Injury = hm_meta$injury[match(colnames(hm_mat), hm_meta$sample_id)],
  Fraction = hm_meta$fraction[match(colnames(hm_mat), hm_meta$sample_id)],
  col = list(Age = c("Young"="#00A087","Aged"="#E64B35"),
             Injury = c("Intact"="#4DBBD5","Fracture"="#F39B7F"),
             Fraction = c("CD45neg"="#B09C85","CD45pos"="#8491B4")))
save_fig("FigS7_DEG_heatmap", 8, 9, {
print(Heatmap(hm_mat, top_annotation = ha, name = "z-score",
              row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
              column_title = "Top 25 up / 25 down aging DEGs"))
})
log_msg("衰老DEG数(", DE_MODE, "): ", nrow(deg_age))
chk <- c("Postn","Dpt","Ctsk","Itm2a","Col1a1","Pdgfra","Csf1r","Il6","Tnf",
         "Spp1","Ccl2","Cxcl12","Cdkn2a","Trp53","Mki67")
print(res_age_df %>% filter(gene %in% chk) %>% select(gene, log2FoldChange, padj))

ego_up <- enrichGO(deg_age %>% filter(log2FoldChange > 0) %>% pull(gene),
                   OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
ego_dn <- enrichGO(deg_age %>% filter(log2FoldChange < 0) %>% pull(gene),
                   OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
write.csv(ego_up@result, "output/tables/TableS2_GO_up.csv", row.names = FALSE)
write.csv(ego_dn@result, "output/tables/TableS3_GO_down.csv", row.names = FALSE)

## ---------------- 7. SenMayo ----------------
senmayo_file <- "data/41467_2022_32552_MOESM4_ESM.xlsx"
if (file.exists(senmayo_file)) {
  senmayo_hs <- readxl::read_xlsx(senmayo_file)[[1]]
} else {
  senmayo_hs <- c("Ccl2","Il6","Serpine1","Icam1","Timp1","Mmp10","Cxcl8","Serpine2",
                  "Tnfsf10c","Plau","Il1b","Timp2","Igfbp7","Mmp1a","Mmp3","Tnf","Fas",
                  "Ccl20","Ccl3","Ccl5","Cxcl1","Cxcl2","Spp1","Angptl4","Gdf15",
                  "Igfbp3","Igfbp5","Igfbp6","Pappa","Stc1","Mif","Hmgb1","Hmgb2",
                  "Anxa1","Vamp3","Vamp5","Vamp7","Arhgdib","Capg","Cops5")
  log_msg("!! 未找到SenMayo补充表，退回核心基因集")
}
map_to_human <- function(ms_genes) {
  og <- safe_orthologs(ms_genes, species = "mouse")
  if (!is.null(og) && nrow(og)) return(unique(na.omit(og$human_symbol)))
  toupper(ms_genes)
}
map_to_mouse <- function(hs_genes) {
  paste0(toupper(substr(hs_genes, 1, 1)), tolower(substr(hs_genes, 2, nchar(hs_genes))))
}
senmayo_ms <- intersect(map_to_mouse(senmayo_hs), rownames(peri))
log_msg("SenMayo小鼠映射: ", length(senmayo_ms), "/", length(senmayo_hs), " 基因可用")
t2g <- data.frame(term = "SenMayo", gene = senmayo_ms)

universe_sc <- res_age_df$gene
ego_sen_up <- enricher(deg_age %>% filter(log2FoldChange > 0) %>% pull(gene),
                       TERM2GENE = t2g, universe = universe_sc, pAdjustMethod = "BH")
gene_list <- res_age_df$stat; names(gene_list) <- res_age_df$gene
gene_list <- sort(gene_list, decreasing = TRUE)
set.seed(123)
gsea_sen <- tryCatch(
  GSEA(gene_list, TERM2GENE = t2g, pAdjustMethod = "BH", minGSSize = 5,
       pvalueCutoff = 1, eps = 1e-10),
  error = function(e) { log_msg("GSEA(SenMayo) 失败: ", e$message); NULL })
gsea_out <- if (is.null(gsea_sen)) data.frame() else gsea_sen@result
for (cc in c("pvalue","p.adjust","qvalue"))
  if (cc %in% names(gsea_out) && is.numeric(gsea_out[[cc]]))
    gsea_out[[cc]] <- formatC(gsea_out[[cc]], format = "e", digits = 3)
write.csv(gsea_out, "output/tables/TableS4_GSEA_SenMayo.csv", row.names = FALSE)
if (gsea_ok(gsea_sen)) {
  log_msg("GSEA(SenMayo): NES=", round(gsea_sen@result$NES[1], 2),
          " p.adj=", signif(gsea_sen@result$p.adjust[1], 3))
  save_fig("Fig2b_GSEA_SenMayo_curve", 7, 5.5, {
  print(gseaplot2(gsea_sen, geneSetID = 1,
                  title = paste0("SenMayo | NES = ", round(gsea_sen@result$NES[1], 2),
                                 ", p.adj = ", formatC(gsea_sen@result$p.adjust[1], format = "e", digits = 2))))
  })
} else {
  log_msg("!! GSEA(SenMayo) 无可用结果（见 TableS4）")
}

peri <- AddModuleScore(peri, features = list(senmayo_ms), name = "SenMayo")
sen_by_donor <- peri@meta.data %>%
  filter(fraction == "CD45neg") %>%
  group_by(donor_id, age_group, injury) %>%
  summarise(SenMayo = mean(SenMayo1), .groups = "drop") %>%
  mutate(age_group = factor(age_group, levels = c("Young","Aged")))
sen_intact <- sen_by_donor %>% filter(injury == "Intact")
sen_test <- tryCatch(wilcox.test(SenMayo ~ age_group, data = sen_intact), error = function(e) NULL)
med_diff <- median(sen_intact$SenMayo[sen_intact$age_group == "Aged"]) /
            median(sen_intact$SenMayo[sen_intact$age_group == "Young"])
auc_val <- if (all(table(sen_intact$age_group) >= 2))
  tryCatch(as.numeric(pROC::auc(sen_intact$age_group, sen_intact$SenMayo)), error = function(e) NA_real_) else NA_real_
sen_cell <- peri@meta.data %>% filter(fraction == "CD45neg", injury == "Intact",
                                      celltype %in% STROMAL)
sen_cell_test <- wilcox.test(SenMayo1 ~ age_group, data = sen_cell)
sen_dir <- sen_cell %>% group_by(age_group) %>%
  summarise(mean_sen = mean(SenMayo1), median_sen = median(SenMayo1), n_cells = n(), .groups = "drop")
print(as.data.frame(sen_dir))
log_msg("SenMayo细胞级方向: Young均值 = ", round(sen_dir$mean_sen[sen_dir$age_group=="Young"], 4),
        " vs Aged均值 = ", round(sen_dir$mean_sen[sen_dir$age_group=="Aged"], 4),
        " → ", ifelse(sen_dir$median_sen[sen_dir$age_group=="Aged"] >
                      sen_dir$median_sen[sen_dir$age_group=="Young"], "Aged>Young(方向符合预期)", "Aged≤Young(需复核)"))
log_msg("SenMayo(CD45neg完整): 供体级中位倍数 = ", round(med_diff, 2),
        if (!is.null(sen_test)) paste0(" | 供体级wilcox p = ", signif(sen_test$p.value, 3)) else " | 供体级n=1v1",
        if (!is.na(auc_val)) paste0(" | AUC(供体级) = ", round(auc_val, 2)) else " | AUC不报告(n=1v1)",
        " | 细胞级wilcox p = ", signif(sen_cell_test$p.value, 3), "（细胞数大，探索性）")
auc_txt <- if (!is.na(auc_val)) {
  paste0("AUC = ", round(auc_val, 2), ", fold = ", round(med_diff, 2), "×")
} else {
  paste0("donor n = 1 vs 1; fold = ", round(med_diff, 2), "× (descriptive)")
}
sen_cell_st <- sen_cell %>% mutate(age_group = factor(age_group, levels = c("Young","Aged")))
p_sen_cell <- ggplot(sen_cell_st, aes(age_group, SenMayo1, fill = age_group)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.1, fill = "white") +
  scale_fill_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
  theme(legend.position = "none") +
  labs(title = "Cell-level SenMayo", x = NULL, y = "SenMayo score",
       subtitle = paste0("CD45neg intact, stromal (n = ", nrow(sen_cell_st), " cells)"))
p_sen_donor <- mk_donor_panel(sen_intact, ylab = "Mean SenMayo score",
                              title = "Donor-level SenMayo") +
  labs(subtitle = auc_txt)
save_fig("Fig2_senmayo_score", 11, 5, {
print(wrap_plots(p_sen_cell, p_sen_donor, ncol = 2))
})
save_fig("FigS2_senmayo_by_celltype", 8, 4.5, {
print(peri@meta.data %>% filter(fraction == "CD45neg") %>%
  mutate(age_group = factor(age_group, levels = c("Young","Aged"))) %>%
  ggplot(aes(reorder(celltype, SenMayo1, median), SenMayo1, fill = age_group)) +
  geom_boxplot(outlier.size = 0.1) + coord_flip() +
  scale_fill_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
  labs(x = NULL, y = "SenMayo score", fill = NULL) +
  theme(plot.margin = margin(r = 45)))
})
write.csv(sen_by_donor, "output/tables/TableS5_SenMayo_by_donor.csv", row.names = FALSE)

## ---------------- 8. age × injury 交互 + 再生响应程序（故事线②） ----------------
peri_neg <- subset(peri, subset = fraction == "CD45neg" & celltype %in% STROMAL)
n_by_grp <- peri_neg@meta.data %>% distinct(donor_id, age_group, injury) %>%
  count(age_group, injury) %>% pivot_wider(names_from = injury, values_from = n, values_fill = 0)
log_msg("CD45neg 每组【文库】数: "); print(as.data.frame(n_by_grp))

if (all(n_by_grp %>% select(-age_group) >= 2)) {
  pb_all <- AggregateExpression(peri_neg, group.by = "sample_id", return.seurat = TRUE)
  mt_all <- peri_neg@meta.data %>% distinct(sample_id, .keep_all = TRUE)
  pb_all <- set_pb_meta(pb_all, mt_all, by = "sample_id")
  pb_all$age_group <- factor(pb_all$age_group, levels = c("Young","Aged"))
  pb_all$injury <- factor(pb_all$injury, levels = c("Intact","Fracture"))
  dds_ix <- run_step("DESeq2_ageXinjury", {
    DESeqDataSetFromMatrix(GetAssayData(pb_all, layer = "counts"),
                           pb_all@meta.data, design = ~ age_group * injury) %>% DESeq()
  })
  res_ix <- results(dds_ix, name = "age_groupAged.injuryFracture")
  write.csv(as.data.frame(res_ix) %>% rownames_to_column("gene"),
            "output/tables/Table2_ageXinjury_interaction.csv", row.names = FALSE)
  log_msg("age×injury 交互: FDR<0.05 = ", sum(res_ix$padj < 0.05, na.rm = TRUE),
          "；名义p<0.01 = ", sum(res_ix$pvalue < 0.01, na.rm = TRUE))
  high_conf <- inner_join(
    deg_age %>% select(gene, sc_lfc = log2FoldChange),
    as.data.frame(res_ix) %>% rownames_to_column("gene") %>%
      filter(pvalue < 0.01) %>% select(gene, slope_inter = log2FoldChange, p_inter = pvalue),
    by = "gene") %>% filter(sign(sc_lfc) == sign(slope_inter))
} else {
  log_msg("!! 每组<2文库 → 跳过样本级交互；high_conf改用细胞级Top50（探索性）")
  res_ix <- NULL
  high_conf <- deg_age %>% filter(abs(log2FoldChange) > 0.25) %>%
    arrange(padj) %>% head(50) %>%
    transmute(gene, sc_lfc = log2FoldChange, slope_inter = NA_real_, p_inter = padj)
}
write.csv(high_conf, "output/tables/Table3_high_confidence_genes.csv", row.names = FALSE)
log_msg("高置信基因数: ", nrow(high_conf))

resp_y <- FindMarkers(subset(peri_neg, age_group == "Young"),
                      ident.1 = "Fracture", ident.2 = "Intact", group.by = "injury",
                      logfc.threshold = 0, min.pct = 0.01, test.use = "wilcox") %>%
  rownames_to_column("gene")
resp_a <- FindMarkers(subset(peri_neg, age_group == "Aged"),
                      ident.1 = "Fracture", ident.2 = "Intact", group.by = "injury",
                      logfc.threshold = 0, min.pct = 0.01, test.use = "wilcox") %>%
  rownames_to_column("gene")
top_y <- resp_y %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  arrange(desc(avg_log2FC)) %>% head(200) %>% pull(gene)
t2g_resp <- data.frame(term = "YoungInjuryResponse", gene = top_y)
gl_a <- resp_a %>% mutate(stat = sign(avg_log2FC) * -log10(p_val + 1e-300)) %>%
  { x <- .$stat; names(x) <- .$gene; x } %>% sort(decreasing = TRUE)
set.seed(123)
gsea_resp <- tryCatch(
  GSEA(gl_a, TERM2GENE = t2g_resp, pAdjustMethod = "BH", minGSSize = 10,
       pvalueCutoff = 1, eps = 1e-10),
  error = function(e) { log_msg("GSEA(youngResponse) 失败: ", e$message); NULL })
write.csv(if (is.null(gsea_resp)) data.frame() else gsea_resp@result,
          "output/tables/TableS6_GSEA_youngResponse_inAged.csv", row.names = FALSE)
if (gsea_ok(gsea_resp)) {
  log_msg("青年损伤响应程序(Top", length(top_y), ")在老年骨折骨膜中GSEA: NES=",
          round(gsea_resp@result$NES[1], 2), " p.adj=", signif(gsea_resp@result$p.adjust[1], 3))
} else {
  log_msg("!! 青年损伤响应程序 GSEA 无可用结果（见 TableS6）")
}
p_fr_cell <- peri@meta.data %>%
  filter(fraction == "CD45neg", injury == "Fracture", celltype %in% STROMAL) %>%
  mutate(age_group = factor(age_group, levels = c("Young","Aged"))) %>%
  ggplot(aes(age_group, SenMayo1, fill = age_group)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.1, fill = "white") +
  scale_fill_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
  theme(legend.position = "none") +
  labs(title = "Cell-level SenMayo after fracture", x = NULL, y = "SenMayo score")
p_fr_donor <- mk_donor_panel(sen_by_donor %>% filter(injury == "Fracture"),
                             ylab = "Mean SenMayo score", title = "Donor-level (n=1v1)")
save_fig("Fig3a_senmayo_after_fracture", 11, 5, {
print(wrap_plots(p_fr_cell, p_fr_donor, ncol = 2))
})
save_fig("Fig3b_gsea_youngResponse", 7, 5, {
  if (gsea_ok(gsea_resp)) {
    print(gseaplot2(gsea_resp, geneSetID = 1,
                    title = paste0("Young injury response in aged | NES = ",
                                   round(gsea_resp@result$NES[1], 2))))
  } else {
    plot.new(); text(0.5, 0.5, "GSEA not available")
  }
})
amp <- inner_join(resp_y %>% select(gene, lfc_y = avg_log2FC),
                  resp_a %>% select(gene, lfc_a = avg_log2FC), by = "gene") %>%
  filter(gene %in% top_y)
amp_test <- wilcox.test(amp$lfc_y, amp$lfc_a, paired = TRUE)
log_msg("青年响应基因诱导幅度: Young中位 = ", round(median(amp$lfc_y), 2),
        " vs Aged中位 = ", round(median(amp$lfc_a), 2),
        " | 配对wilcox p = ", signif(amp_test$p.value, 3))

## 8b. GSE198666：骨折骨痂 年轻 vs 老龄（同物种验证）
run_step("GSE198666_callus", {
  if (!DO_CALLUS_VALIDATION) { log_msg("跳过"); return(NULL) }
  fetch_198 <- fetch_geo("GSE198666", "data/GSE198666")
  pd_198 <- fetch_198$pd
  meta_198 <- tibble(gsm = pd_198$geo_accession, title = pd_198$title,
                     source = pd_198$source_name_ch1,
                     txt = tolower(paste(title, source))) %>%
    mutate(age_group = case_when(
      str_detect(txt, "young|\\b[2-6] ?mo|adult") ~ "Young",
      str_detect(txt, "aged|old|\\b(1[89]|2[0-9]) ?mo") ~ "Aged",
      TRUE ~ NA_character_))
  if (any(is.na(meta_198$age_group))) stop("GSE198666 分组解析失败 → 请override")
  meta_198 <- meta_198 %>%
    mutate(donor_id = gsm, sample_id = gsm, injury = "Fracture", fraction = "CD45pos")
  objs198 <- setNames(lapply(seq_len(nrow(meta_198)), function(i)
    read_sample_any(fetch_198$rawdir, meta_198$gsm[i], meta_198[i, ], mt.pattern = "^mt-")),
    meta_198$gsm)
  seu198 <- reduce(objs198, merge) %>% NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = "donor_id") %>%
    FindNeighbors(reduction = "harmony") %>% FindClusters(resolution = 0.5) %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
  ## v3.8【补丁4】strengthening 脚本 §7 时间序列输入
  write.csv(meta_198, "output/tables/TableS0_meta198.csv", row.names = FALSE)
  saveRDS(seu198, "output/seu198.rds")
  log_msg("补丁4: TableS0_meta198.csv / seu198.rds 已保存")
  pb198 <- AggregateExpression(seu198, group.by = "donor_id", return.seurat = TRUE)
  mt198 <- seu198@meta.data %>% distinct(donor_id, .keep_all = TRUE)
  pb198 <- set_pb_meta(pb198, mt198, by = "donor_id")
  pb198$age_group <- factor(pb198$age_group, levels = c("Young","Aged"))
  dds198 <- DESeqDataSetFromMatrix(GetAssayData(pb198, layer = "counts"),
                                   pb198@meta.data, design = ~ age_group) %>% DESeq()
  res198 <- results(dds198, contrast = c("age_group","Aged","Young"))
  res198_df <- as.data.frame(res198) %>% rownames_to_column("gene")
  write.csv(res198_df, "output/tables/Table4b_callus_aging_DEGs.csv", row.names = FALSE)
  cv198 <- inner_join(deg_age %>% select(gene, lfc_peri = log2FoldChange),
                      res198_df %>% select(gene, lfc_callus = log2FoldChange), by = "gene") %>%
    mutate(direction_consistent = sign(lfc_peri) == sign(lfc_callus))
  rho198 <- spearman_rho(cv198$lfc_peri, cv198$lfc_callus)
  consistency_pct <- round(mean(cv198$direction_consistent, na.rm = TRUE) * 100, 1)
  log_msg("骨痂验证: 共有基因 ", nrow(cv198), "，方向一致率 ",
          consistency_pct, "%，rho = ", round(rho198, 2))
  write.csv(cv198, "output/tables/Table4_callus_validation.csv", row.names = FALSE)
  up200 <- deg_age %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% head(200) %>% pull(gene)
  t2g_up <- data.frame(term = "PeriAgingUP", gene = intersect(up200, rownames(res198)))
  gl198 <- res198$stat; names(gl198) <- rownames(res198); gl198 <- sort(gl198, decreasing = TRUE)
  set.seed(123)
  gsea198 <- tryCatch(
    GSEA(gl198, TERM2GENE = t2g_up, pAdjustMethod = "BH", minGSSize = 10,
         pvalueCutoff = 1, eps = 1e-10),
    error = function(e) { log_msg("GSEA(callus) 失败: ", e$message); NULL })
  write.csv(if (is.null(gsea198)) data.frame() else gsea198@result,
            "output/tables/TableS6b_GSEA_periAging_inCallus.csv", row.names = FALSE)
  if (gsea_ok(gsea198)) {
    log_msg("骨膜衰老程序在骨痂中GSEA: NES=", round(gsea198@result$NES[1],2),
            " p.adj=", signif(gsea198@result$p.adjust[1],3))
  } else {
    log_msg("!! 骨膜衰老程序在骨痂中 GSEA 无可用结果（见 TableS6b）")
  }
  save_fig("Fig3c_callus_scatter", 6, 5, {
  print(ggplot(cv198, aes(lfc_peri, lfc_callus, color = direction_consistent)) +
        geom_point(size = 1.5, alpha = .5) +
        geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
        scale_color_manual(values = c("TRUE"="#E64B35","FALSE"="grey70")) +
        labs(x = "Periosteum aging log2FC", y = "Callus aging log2FC", color = NULL) +
        ggtitle(paste0("Periosteum-callus replication (n=", nrow(cv198), ")"),
                subtitle = paste0("Direction consistency = ", consistency_pct, "%")) +
        theme(legend.position = "none"))
  })
  save_fig("Fig3d_callus_gsea", 7, 5, {
    if (gsea_ok(gsea198)) {
      print(gseaplot2(gsea198, geneSetID = 1,
                      title = paste0("Peri aging-UP in callus | NES = ",
                                     round(gsea198@result$NES[1], 2))))
    } else {
      plot.new(); text(0.5, 0.5, "GSEA not available")
    }
  })
})

## ---------------- 9. 独立图谱验证（GSE297256）+ 人类验证层 ----------------
## v3.8【重写】主路径：解包目录名硬解析 Y/M/O + 直接 ReadMtx + 免 harmony
##      （修复 v3.7 实测 "invalid 'data'"）；title 解析 + read_sample_any 保留为回退
run_step("GSE297256_atlas", {
  if (!DO_ATLAS_VALIDATION) { log_msg("跳过"); return(NULL) }
  fetch_297 <- fetch_geo("GSE297256", "data/GSE297256")
  pd_297 <- fetch_297$pd

  untar_dirs <- list.dirs(fetch_297$rawdir, recursive = TRUE) %>%
    keep(~ str_detect(basename(.x), "^[YMO]_P$"))
  if (length(untar_dirs) >= 3) {
    meta_297 <- tibble(
      folder   = untar_dirs,
      gsm      = basename(dirname(untar_dirs)) %>% str_extract("GSM[0-9]+"),
      tag      = substr(basename(untar_dirs), 1, 1),
      age_mo   = recode(tag, Y = 3, M = 9, O = 18),
      donor_id = gsm, sample_id = gsm,
      age_group = "Mixed", injury = "Intact", fraction = "Mixed")
    if (any(is.na(meta_297$gsm)) || any(is.na(meta_297$age_mo))) {
      print(as.data.frame(meta_297)); stop("v3.8 目录名解析失败 → 请检查 data/GSE297256/RAW 结构") }
    log_msg("v3.8: GSE297256 元数据来自解包目录名 (Y/M/O)")
    print(as.data.frame(meta_297 %>% select(gsm, tag, age_mo)))
    objs297 <- setNames(lapply(seq_len(nrow(meta_297)), function(i) {
      m <- meta_297[i, ]
      counts <- ReadMtx(mtx      = file.path(m$folder, "matrix.mtx.gz"),
                        cells    = file.path(m$folder, "barcodes.tsv.gz"),
                        features = file.path(m$folder, "features.tsv.gz"),
                        feature.column = 2)
      obj <- CreateSeuratObject(counts, project = m$sample_id,
                                min.cells = 3, min.features = 200)
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
      obj$age_mo    <- m$age_mo
      obj$donor_id  <- m$donor_id
      obj$sample_id <- m$sample_id
      obj$age_group <- "Mixed"; obj$injury <- "Intact"; obj$fraction <- "Mixed"
      subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
    }), meta_297$sample_id)
    ## 3 文库 × 3 donor：样本量小，直接 merge 降维即可（免 harmony，更稳）
    seu297 <- reduce(objs297, merge) %>% NormalizeData() %>% FindVariableFeatures() %>%
      ScaleData() %>% RunPCA() %>%
      FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.5) %>%
      RunUMAP(dims = 1:30)
  } else {
    ## ---- 回退路径：v3.7b title/source 解析 + read_sample_any ----
    age_tag_from_raw <- function(rawdir, gsm) {
      ff <- list.files(rawdir, pattern = gsm, recursive = TRUE)
      if (!length(ff)) return(NA_character_)
      tag <- str_match(basename(ff[1]), "_([YMO])[_-]P")[, 2]
      tag[1]
    }
    meta_297 <- tibble(
      gsm = pd_297$geo_accession, title = pd_297$title,
      source = pd_297$source_name_ch1,
      txt = tolower(paste(title, source))) %>%
      mutate(
        age_mo = suppressWarnings(as.numeric(str_extract(
          txt, "(?<![0-9])(3|9|18) ?(mo|m)(?= |$|[^0-9])"))),
        age_mo = case_when(!is.na(age_mo) ~ age_mo,
          str_detect(txt, "young|(^|[^a-z])y[ _-]?p([^a-z]|$)") ~ 3,
          str_detect(txt, "mid|middle|(^|[^a-z])m[ _-]?p([^a-z]|$)") ~ 9,
          str_detect(txt, "old|aged|(^|[^a-z])o[ _-]?p([^a-z]|$)") ~ 18,
          TRUE ~ NA_real_),
        age_mo = case_when(!is.na(age_mo) ~ age_mo,
          age_tag_from_raw(fetch_297$rawdir, gsm) == "Y" ~ 3,
          age_tag_from_raw(fetch_297$rawdir, gsm) == "M" ~ 9,
          age_tag_from_raw(fetch_297$rawdir, gsm) == "O" ~ 18,
          TRUE ~ NA_real_),
        age_group = case_when(age_mo <= 4 ~ "Young", age_mo <= 12 ~ "Mid", TRUE ~ "Aged"),
        donor_id = gsm, sample_id = gsm, injury = "Intact", fraction = "Mixed")
    if (any(is.na(meta_297$age_mo))) {
      print(as.data.frame(meta_297)); stop("GSE297256 年龄解析失败 → 请override") }
    objs297 <- setNames(lapply(seq_len(nrow(meta_297)), function(i)
      read_sample_any(fetch_297$rawdir, meta_297$gsm[i], meta_297[i, ], mt.pattern = "^mt-")),
      meta_297$gsm)
    seu297 <- reduce(objs297, merge) %>% NormalizeData() %>% FindVariableFeatures() %>%
      ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = "donor_id") %>%
      FindNeighbors(reduction = "harmony") %>% FindClusters(resolution = 0.5) %>%
      RunUMAP(reduction = "harmony", dims = 1:30)
  }
  pb297 <- AggregateExpression(seu297, group.by = "donor_id", return.seurat = TRUE)
  mt297 <- seu297@meta.data %>% distinct(donor_id, .keep_all = TRUE)
  pb297 <- set_pb_meta(pb297, mt297, by = "donor_id")
  pb297$age_mo <- as.numeric(pb297$age_mo)
  stopifnot(all(!is.na(pb297$age_mo)))
  dds297 <- DESeqDataSetFromMatrix(GetAssayData(pb297, layer = "counts"),
                                   pb297@meta.data, design = ~ age_mo) %>% DESeq()
  res297 <- results(dds297, name = "age_mo")
  res297_df <- as.data.frame(res297) %>% rownames_to_column("gene")
  write.csv(res297_df, "output/tables/Table4c_atlas_age_slope.csv", row.names = FALSE)
  ## v3.8【补丁3】strengthening 脚本 §8 衰老时钟训练集
  save(pb297, file = "output/pb297.RData")
  log_msg("补丁3: pb297.RData 已保存")
  cv297 <- inner_join(deg_age %>% select(gene, lfc_peri = log2FoldChange),
                      res297_df %>% select(gene, slope_atlas = log2FoldChange), by = "gene") %>%
    mutate(direction_consistent = sign(lfc_peri) == sign(slope_atlas))
  rho297 <- spearman_rho(cv297$lfc_peri, cv297$slope_atlas)
  log_msg("GSE297256独立验证: 共有基因 ", nrow(cv297), "，方向一致率 ",
          round(mean(cv297$direction_consistent, na.rm = TRUE)*100,1), "%，rho = ", round(rho297,2))
  write.csv(cv297, "output/tables/Table4d_atlas_validation.csv", row.names = FALSE)
  if (nrow(cv297) > 10) {
    save_fig("Fig4_atlas_validation", 6, 5.5, {
    print(ggplot(cv297, aes(lfc_peri, slope_atlas, color = direction_consistent)) +
          geom_point(size = 1.5, alpha = .5) +
          geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
          scale_color_manual(values = c("TRUE"="#E64B35","FALSE"="grey70")) +
          labs(x = "GSE280914 aging log2FC (Aged vs Young)",
               y = "GSE297256 age slope (per month)",
               title = paste0("Independent atlas replication (n=", nrow(cv297), ")"),
               subtitle = paste0("consistency = ",
                                 round(mean(cv297$direction_consistent, na.rm=TRUE)*100,1), "%")) +
          theme(legend.position = "none"))
    })
  }
})

run_step("GSE232516", {
  fetch_232 <- fetch_geo("GSE232516", "data/GSE232516")
  pd_232 <- fetch_232$pd
  meta_232 <- tibble(gsm = pd_232$geo_accession, title = pd_232$title,
                     source = pd_232$source_name_ch1,
                     txt = tolower(paste(title, source))) %>%
    mutate(group = case_when(
      str_detect(txt, "cpt|pseudarth|patholog|tibia") ~ "CPT",
      str_detect(txt, "iliac|crest|control|normal|healthy") ~ "Control",
      TRUE ~ NA_character_))
  if (any(is.na(meta_232$group))) stop("GSE232516 分组解析失败 → 请override")
  meta_232 <- meta_232 %>%
    mutate(donor_id = gsm, sample_id = gsm, age_group = group, injury = "Intact", fraction = "Mixed")
  objs232 <- setNames(lapply(seq_len(nrow(meta_232)), function(i)
    read_sample_any(fetch_232$rawdir, meta_232$gsm[i], meta_232[i, ], mt.pattern = "^MT-") ),
    meta_232$gsm)
  seu232 <- reduce(objs232, merge) %>% NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = "donor_id") %>%
    FindNeighbors(reduction = "harmony") %>% FindClusters(resolution = 0.5) %>%
    RunUMAP(reduction = "harmony", dims = 1:30)
  hum_markers <- c("PDGFRA","COL1A1","PRRX1","CTSK","POSTN","DPT","ITM2A","LRP1","ANPEP",
                   "PTPRC","CD3E","CSF1R","CD14","PECAM1","KDR","RGS5","MYOG")
  hum_markers <- hum_markers[hum_markers %in% rownames(seu232)]
  save_fig("FigS3_232_dotplot", 10, 6, { print(DotPlot(seu232, features = hum_markers) + RotatedAxis()) })
  seu232$celltype <- case_when(seu232$seurat_clusters %in% c("0","1") ~ "pSSPC", TRUE ~ "Other")
  pb232 <- AggregateExpression(seu232, group.by = "donor_id", return.seurat = TRUE)
  mt232 <- seu232@meta.data %>% distinct(donor_id, .keep_all = TRUE)
  pb232 <- set_pb_meta(pb232, mt232, by = "donor_id")
  pb232$group <- factor(pb232$age_group, levels = c("Control","CPT"))
  stopifnot(all(!is.na(pb232$group)))
  dds232 <- DESeqDataSetFromMatrix(GetAssayData(pb232, layer = "counts"),
                                   pb232@meta.data, design = ~ group) %>% DESeq()
  res232 <- results(dds232, contrast = c("group","CPT","Control"))
  ## v3.8: babelgene 失败时回退 toupper 直映射（探索性可接受）
  og_map <- safe_orthologs(deg_age$gene, species = "mouse")
  if (!is.null(og_map)) {
    deg_age_hs <- deg_age %>% inner_join(og_map, by = c("gene" = "symbol")) %>%
      select(gene_hs = human_symbol, lfc_ms = log2FoldChange)
  } else {
    deg_age_hs <- deg_age %>% transmute(gene_hs = toupper(gene), lfc_ms = log2FoldChange)
    log_msg("v3.8: 跨物种映射使用 toupper 直映射（方法学限制中需声明）")
  }
  res232_df <- as.data.frame(res232) %>% rownames_to_column("gene_hs")
  cross_sp <- inner_join(deg_age_hs, res232_df, by = "gene_hs") %>%
    mutate(direction_consistent = sign(lfc_ms) == sign(log2FoldChange))
  rho_sp <- spearman_rho(cross_sp$lfc_ms, cross_sp$log2FoldChange)
  write.csv(cross_sp, "output/tables/Table5_cross_species_validation.csv", row.names = FALSE)
  log_msg("人骨膜验证: 共有基因 ", nrow(cross_sp), "，方向一致率 ",
          round(mean(cross_sp$direction_consistent, na.rm=TRUE)*100, 1), "%，rho = ", round(rho_sp, 2))
  if (nrow(cross_sp) > 10) {
    save_fig("Fig5_cross_species_scatter", 6, 5.5, {
    print(ggplot(cross_sp, aes(lfc_ms, log2FoldChange, color = direction_consistent)) +
          geom_point(size = 2, alpha = .7) +
          geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
          geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
          scale_color_manual(values = c("TRUE"="#E64B35","FALSE"="grey70"),
                             labels = c("TRUE"="Consistent","FALSE"="Discordant")) +
          labs(x = "Mouse aging log2FC (Aged vs Young)",
               y = "Human CPT vs Control log2FC",
               title = paste0("Cross-species replication (n = ", nrow(cross_sp), ")"),
               subtitle = paste0("direction consistency = ",
                                 round(mean(cross_sp$direction_consistent, na.rm=TRUE)*100,1), "%")))
    })
  }
})

## ---------------- 9b. GSE278165 跨物种锚定（人 P-SSC） ----------------
run_step("GSE278165", {
  if (!DO_HUMAN_ATLAS) { log_msg("DO_HUMAN_ATLAS=FALSE，跳过"); return(NULL) }
  fetch_278 <- fetch_geo("GSE278165", "data/GSE278165")
  pd_278 <- fetch_278$pd
  meta_278 <- tibble(gsm = pd_278$geo_accession, title = pd_278$title,
                     source = pd_278$source_name_ch1) %>%
    mutate(age_group = "Human", injury = "Intact", fraction = "Mixed",
           donor_id = gsm, sample_id = gsm)
  rds_files <- list.files(fetch_278$rawdir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  if (length(rds_files)) {
    hum <- readRDS(rds_files[1])
  } else {
    objs278 <- setNames(lapply(seq_len(nrow(meta_278)), function(i)
      read_sample_any(fetch_278$rawdir, meta_278$gsm[i], meta_278[i, ], mt.pattern = "^MT-")),
      meta_278$gsm)
    hum <- reduce(objs278, merge) %>% NormalizeData() %>% FindVariableFeatures()
  }
  if (ncol(hum) > 20000) hum <- subset(hum, cells = sample(colnames(hum), 20000))
  if (!"seurat_clusters" %in% colnames(hum@meta.data)) {
    hum <- hum %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
      FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.5)
  }
  ## v3.8: babelgene 失败时回退 toupper
  og2 <- safe_orthologs(rownames(peri), species = "mouse")
  if (is.null(og2)) {
    log_msg("v3.8: 无可用同源映射 → 跳过 GSE278165 锚定（不影响主线）"); return(NULL) }
  ms2hs <- og2 %>% filter(!is.na(human_symbol)) %>%
    group_by(symbol) %>% slice(1) %>% ungroup()
  counts_ms <- GetAssayData(peri, assay = "RNA", layer = "counts")
  shared <- intersect(rownames(counts_ms), ms2hs$symbol)
  mat_hs <- counts_ms[shared, , drop = FALSE]
  rownames(mat_hs) <- ms2hs$human_symbol[match(shared, ms2hs$symbol)]
  dup <- duplicated(rownames(mat_hs))
  mat_hs <- mat_hs[!dup, , drop = FALSE]
  ms_hs <- CreateSeuratObject(mat_hs, meta.data = peri@meta.data[, c("age_group","injury","celltype")])
  ms_hs <- NormalizeData(ms_hs) %>% FindVariableFeatures()
  hum <- FindVariableFeatures(NormalizeData(hum))
  anchors <- FindTransferAnchors(reference = hum, query = ms_hs, dims = 1:30, reduction = "cca")
  td <- TransferData(anchorset = anchors, refdata = hum$seurat_clusters,
                     dims = 1:30, weight.reduction = "cca")
  ms_hs <- AddMetaData(ms_hs, td)
  ms_hs$predicted.hum_label <- ms_hs$predicted.id
  pred_df <- ms_hs@meta.data %>%
    mutate(age_group = factor(age_group, levels = c("Young","Aged"))) %>%
    count(age_group, injury, celltype, predicted.hum_label) %>%
    group_by(age_group, injury, celltype) %>% mutate(pct = n/sum(n)*100)
  write.csv(pred_df, "output/tables/TableS7_label_transfer.csv", row.names = FALSE)
  save_fig("FigS4_label_transfer", 8, 5, {
  print(ggplot(pred_df, aes(celltype, pct, fill = predicted.hum_label)) +
        geom_col() + coord_flip() + facet_wrap(~ age_group) +
        labs(x = NULL, y = "%", fill = "Human cluster"))
  })
})

## ---------------- 10. CellChat：青年 vs 老龄 ----------------
library(CellChat); library(tidyverse); library(Matrix)
data(CellChatDB.mouse)

peri_cc <- subset(peri, subset = injury == "Intact")
meta_cc <- peri_cc@meta.data %>% rownames_to_column("cell") %>% mutate(labels = celltype)
set.seed(123)
cells_cc <- meta_cc %>% group_by(age_group, labels) %>% slice_sample(n = 4000) %>% pull(cell)
cc_labels  <- meta_cc$labels[match(cells_cc, meta_cc$cell)]
cc_age     <- meta_cc$age_group[match(cells_cc, meta_cc$cell)]
cc_raw <- subset(peri_cc, cells = cells_cc)
ly <- grep("^data", Layers(cc_raw, assay = "RNA"), value = TRUE)
mats <- lapply(ly, function(l) GetAssayData(cc_raw, assay = "RNA", layer = l))
common_genes <- Reduce(intersect, lapply(mats, rownames))
data_in <- do.call(cbind, lapply(mats, function(m) m[common_genes, , drop = FALSE]))
rm(mats, cc_raw); gc()
sig_genes <- unique(unlist(strsplit(c(as.character(CellChatDB.mouse$interaction$ligand),
                                      as.character(CellChatDB.mouse$interaction$receptor)), "_")))
gene.keep <- rownames(data_in)[Matrix::rowSums(data_in > 0) > 50]
data_use  <- data_in[union(gene.keep, intersect(sig_genes, rownames(data_in))), ]
rm(data_in); gc()
ord <- match(colnames(data_use), cells_cc)
stopifnot(!any(is.na(ord)))
meta_cc2 <- data.frame(labels = cc_labels[ord], age = cc_age[ord],
                       row.names = colnames(data_use), stringsAsFactors = FALSE)
rm(list = intersect(c("peri_cc","peri_intact","pb","objs","peri_neg","pb_all"), ls()))
gc(reset = TRUE)

run_cellchat <- function(data_mat, meta, group.by = "labels") {
  cc <- createCellChat(object = as.matrix(data_mat), meta = meta, group.by = group.by)
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc) %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(population.size = TRUE) %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  cc
}
idx_Y <- meta_cc2$age == "Young"; idx_A <- meta_cc2$age == "Aged"
log_msg("CellChat Young ...")
cc_Y <- run_cellchat(data_use[, idx_Y], meta_cc2[idx_Y, ])
log_msg("CellChat Aged ...")
cc_A <- run_cellchat(data_use[, idx_A], meta_cc2[idx_A, ])
save(cc_Y, cc_A, meta_cc2, file = "output/CellChat_periosteum.RData")
object.list <- list(Young = cc_Y, Aged = cc_A)
merged.cc <- mergeCellChat(object.list, add.names = names(object.list))

write.csv(subsetCommunication(cc_Y), "output/tables/TableS8_CC_Young.csv", row.names = FALSE)
write.csv(subsetCommunication(cc_A), "output/tables/TableS9_CC_Aged.csv", row.names = FALSE)

save_fig("Fig6b_CC_chord_SPP1_TGFB", 18, 9, {
  par(mfrow = c(1, 2), mar = c(1, 1, 2, 1), cex.main = 0.9)
  tryCatch(netVisual_chord_gene(cc_Y, signaling = c("SPP1","TGFb"),
                                lab.cex = 0.4, title.name = "Young (intact)"),
           error = function(e) { log_msg("chord Young 失败: ", e$message)
                                 plot.new(); text(0.5, 0.5, "chord failed (Young)") })
  tryCatch(netVisual_chord_gene(cc_A, signaling = c("SPP1","TGFb"),
                                lab.cex = 0.4, title.name = "Aged (intact)"),
           error = function(e) { log_msg("chord Aged 失败: ", e$message)
                                 plot.new(); text(0.5, 0.5, "chord failed (Aged)") })
  par(mfrow = c(1, 1))
})

save_fig("Fig6c_cell_chord_SPP1_TGFB", 18, 9, {
  par(mfrow = c(1, 2), mar = c(1, 1, 2, 1), cex.main = 0.9)
  netVisual_chord_cell(cc_Y, signaling = c("SPP1","TGFb"),
                       lab.cex = 1.1, title.name = "Young (intact)")
  netVisual_chord_cell(cc_A, signaling = c("SPP1","TGFb"),
                       lab.cex = 1.1, title.name = "Aged (intact)")
  par(mfrow = c(1, 1))
})

df_cc <- tibble(
  age = c("Young","Aged"),
  n_interaction = c(nrow(subsetCommunication(cc_Y)), nrow(subsetCommunication(cc_A))),
  strength = c(sum(cc_Y@net$weight), sum(cc_A@net$weight)))
print(df_cc)
save_fig("Fig6_CC_overview", 9.5, 4.5, {
print(ggplot(df_cc %>% pivot_longer(-age, names_to = "measure", values_to = "value") %>%
               mutate(age = factor(age, levels = c("Young","Aged"))),
             aes(age, value, fill = age)) +
      geom_col(width = .6) +
      geom_text(aes(label = round(value, 3)), vjust = -0.4) +
      facet_wrap(~measure, scales = "free_y") +
      scale_fill_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
      theme(legend.position = "none", axis.title = element_blank()))
})
sen_lig <- c("Csf1","Spp1","Il6","Tnf","Mif","Tgfb1","Cxcl12","Ccl2")
sen_comm <- bind_rows(
  subsetCommunication(cc_A) %>% filter(ligand %in% sen_lig) %>% mutate(age = "Aged"),
  subsetCommunication(cc_Y) %>% filter(ligand %in% sen_lig) %>% mutate(age = "Young")) %>%
  select(age, ligand, receptor, source, target, prob, pval) %>% arrange(desc(prob))
write.csv(sen_comm, "output/tables/TableS9b_CC_senescence_ligands.csv", row.names = FALSE)
log_msg("衰老配体通讯Top边: "); print(head(as.data.frame(sen_comm), 15))
sasp_pathways <- intersect(c("CSF","TNF","IL6","OSM","TGFb","SPP1","CCL","CXCL","MIF","BAFF"),
                           unique(c(cc_A@netP$pathways, cc_Y@netP$pathways)))
save_fig("FigS5_CC_sasp_bubble", 12, 9, {
print(netVisual_bubble(merged.cc, signaling = sasp_pathways, comparison = c(1, 2),
                       angle.x = 90, remove.isolate = TRUE, font.size = 5) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
            legend.position = "right"))
})

## ---------------- 11. Monocle3：激活轨迹停滞 + 衰老耦合 ----------------
library(monocle3)
peri_tr <- subset(peri, subset = fraction == "CD45neg" & celltype %in% STROMAL)
set.seed(123)
meta_tr <- peri_tr@meta.data %>% rownames_to_column("cell") %>% mutate(ct = celltype)
cells_tr <- meta_tr %>% group_by(age_group, injury, ct) %>% slice_sample(n = 1500) %>% pull(cell)
saveRDS(subset(peri_tr, cells = cells_tr), "output/tr_obj.rds")
tr_obj <- readRDS("output/tr_obj.rds")
pd <- tr_obj@meta.data
pd$celltype <- as.character(Idents(tr_obj))
pd$group <- factor(paste(pd$age_group, pd$injury, sep = "_"),
                   levels = c("Young_Intact","Young_Fracture","Aged_Intact","Aged_Fracture"))
ly <- grep("^data", Layers(tr_obj, assay = "RNA"), value = TRUE)
mats <- lapply(ly, function(l) GetAssayData(tr_obj, assay = "RNA", layer = l))
common_genes <- Reduce(intersect, lapply(mats, rownames))
expr <- do.call(cbind, lapply(mats, function(m) m[common_genes, , drop = FALSE]))
rm(mats, tr_obj); gc()
keep_cells <- intersect(cells_tr, colnames(expr))
expr <- expr[, keep_cells]; pd <- pd[keep_cells, ]
tr_gene_pool <- c(deg_age$gene[1:3000], senmayo_ms,
                  c("Prrx1","Pdgfra","Itm2a","Dpt","Postn","Col1a1","Acta2","Sp7","Runx2",
                    "Ibsp","Sox9","Acan","Mki67","Cdkn2a","Trp53"))
tr_genes_use <- intersect(unique(tr_gene_pool), rownames(expr))
expr_tr <- expr[tr_genes_use, ]; rm(expr); gc()
cds <- new_cell_data_set(expression_data = expr_tr, cell_metadata = pd,
                         gene_metadata = data.frame(gene_short_name = rownames(expr_tr),
                                                    row.names = rownames(expr_tr)))
set.seed(123)
cds <- preprocess_cds(cds, num_dim = 30) %>%
  reduce_dimension(reduction_method = "UMAP") %>%
  cluster_cells(reduction_method = "UMAP") %>%
  learn_graph(use_partition = FALSE)
ROOT_CELLS <- colnames(cds)[colData(cds)$celltype == "pSSPC_q" &
                            colData(cds)$age_group == "Young" &
                            colData(cds)$injury == "Intact"][1:20]
stopifnot(length(ROOT_CELLS) > 0)
cds <- order_cells(cds, root_cells = ROOT_CELLS)
if (!"pseudotime" %in% names(colData(cds)))
  colData(cds)$pseudotime <- cds@principal_graph_aux$UMAP$pseudotime
if ("SenMayo1" %in% colnames(pd)) {
  colData(cds)$SenMayo1 <- as.numeric(pd$SenMayo1)
} else {
  sen_in <- intersect(senmayo_ms, rownames(cds))
  stopifnot(length(sen_in) > 0)
  colData(cds)$SenMayo1 <- as.numeric(Matrix::colMeans(assay(cds[sen_in, ], "counts")))
}
save(cds, file = "output/monocle3_cds_periosteum.RData")

pt <- data.frame(pseudotime = colData(cds)$pseudotime,
                 group = colData(cds)$group,
                 age_group = colData(cds)$age_group,
                 injury = colData(cds)$injury,
                 celltype = colData(cds)$celltype,
                 senmayo = colData(cds)$SenMayo1)
write.csv(pt, "output/tables/TableS10_pseudotime.csv", row.names = FALSE)
pt_med <- pt %>% group_by(group) %>%
  summarise(median_pt = median(pseudotime), mean_sen = mean(senmayo), n = n(), .groups = "drop")
print(as.data.frame(pt_med)); write.csv(pt_med, "output/tables/TableS10b_pt_by_group.csv", row.names = FALSE)
stagnation <- pt %>% group_by(age_group, injury) %>%
  summarise(m = median(pseudotime), .groups = "drop") %>%
  pivot_wider(names_from = injury, values_from = m) %>%
  mutate(delta = Fracture - Intact)
print(as.data.frame(stagnation))
yf_af <- list(
  YF_vs_YI = wilcox.test(pseudotime ~ injury, data = pt %>% filter(age_group == "Young")),
  AF_vs_AI = wilcox.test(pseudotime ~ injury, data = pt %>% filter(age_group == "Aged")),
  AF_vs_YF = wilcox.test(pseudotime ~ group, data = pt %>%
                           filter(group %in% c("Young_Fracture","Aged_Fracture"))))
log_msg("推进幅度 ΔYoung = ", round(stagnation$delta[stagnation$age_group=="Young"], 2),
        " vs ΔAged = ", round(stagnation$delta[stagnation$age_group=="Aged"], 2),
        " | YFvsYI p=", signif(yf_af$YF_vs_YI$p.value,3),
        " AFvsAI p=", signif(yf_af$AF_vs_AI$p.value,3),
        " AFvsYF p=", signif(yf_af$AF_vs_YF$p.value,3))
pt_cor <- cor.test(pt$pseudotime, pt$senmayo, method = "spearman")
log_msg("轨迹×衰老耦合: rho = ", round(unname(pt_cor$estimate), 2),
        " (p = ", signif(pt_cor$p.value, 3), ")")

rd <- tryCatch(reducedDims(cds)[["UMAP"]], error = function(e) NULL)
if (is.null(rd) && length(reducedDims(cds)) > 0) rd <- reducedDims(cds)[[1]]
rd <- as.data.frame(as.matrix(rd))
stopifnot(ncol(rd) >= 2)
names(rd)[1:2] <- c("UMAP_1", "UMAP_2")
rd$cell <- rownames(rd)
pt_df <- data.frame(cell = colnames(cds),
                    pseudotime = as.numeric(colData(cds)$pseudotime),
                    celltype = as.character(colData(cds)$celltype),
                    stringsAsFactors = FALSE)
miles <- rd %>% inner_join(pt_df, by = "cell") %>%
  group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),
            pt = median(pseudotime, na.rm = TRUE), .groups = "drop") %>%
  arrange(pt) %>% mutate(order = seq_len(n()))
root_xy <- rd %>% filter(cell %in% ROOT_CELLS[ROOT_CELLS %in% rd$cell]) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
milestone_layer <- list()
if (nrow(miles) >= 2) {
  seg <- miles %>% mutate(xend = dplyr::lead(UMAP_1), yend = dplyr::lead(UMAP_2)) %>%
    filter(!is.na(xend))
  milestone_layer <- c(milestone_layer, list(
    geom_segment(data = seg, aes(x = UMAP_1, y = UMAP_2, xend = xend, yend = yend),
                 arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"),
                 color = "grey45", linewidth = 0.4, inherit.aes = FALSE)))
}
if (requireNamespace("ggrepel", quietly = TRUE)) {
  milestone_layer <- c(milestone_layer, list(
    ggrepel::geom_label_repel(data = miles, aes(UMAP_1, UMAP_2, label = order),
                              size = 3, fontface = "bold", fill = "white",
                              label.size = 0.2, box.padding = 0.35, point.padding = 0.25,
                              inherit.aes = FALSE, show.legend = FALSE)))
} else {
  milestone_layer <- c(milestone_layer, list(
    geom_label(data = miles, aes(UMAP_1, UMAP_2, label = order), size = 3,
               fontface = "bold", inherit.aes = FALSE, show.legend = FALSE)))
}
milestone_layer <- c(milestone_layer, list(
  geom_point(data = root_xy, aes(UMAP_1, UMAP_2), shape = 8, size = 4,
             color = "#E64B35", inherit.aes = FALSE, show.legend = FALSE)))

p7_1 <- plot_cells(cds, color_cells_by = "celltype", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5)
p7_2 <- plot_cells(cds, color_cells_by = "group", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5)
p7_3 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5)
for (lyr in milestone_layer) p7_3 <- p7_3 + lyr
p7_3 <- p7_3 + labs(title = "Pseudotime (numbered order)",
                    caption = "★ = root (Young Intact pSSPC_q); 1..n = cell-state milestones by median pseudotime")
save_fig("Fig7_trajectory", 14, 4.8, {
print(wrap_plots(p7_1, p7_2, p7_3, ncol = 3))
})
save_fig("Fig8_pt_density_sencoupling", 11, 4.5, {
print(wrap_plots(
  ggplot(pt, aes(pseudotime, fill = group)) +
    geom_density(alpha = .6) + labs(x = "Pseudotime", y = "Density", fill = NULL),
  ggplot(pt, aes(pseudotime, senmayo, color = group)) +
    geom_point(size = .6, alpha = .4) + geom_smooth(method = "loess", se = FALSE) +
    labs(title = paste0("SenMayo along pseudotime | rho = ",
                        round(unname(pt_cor$estimate), 2)),
         x = "Pseudotime", y = "SenMayo", color = NULL),
  ncol = 2))
})
save(pt, pt_med, stagnation, pt_cor, yf_af, file = "output/pseudotime_results.RData")

## ---------------- 12. cis-MR（故事线③） ----------------
run_step("MR", {
  if (length(OUTCOME_IDS) == 0)
    stop("请先运行 pick_outcomes('fracture') 查询真实结局ID，填入脚本顶部 OUTCOME_IDS")
  library(TwoSampleMR); library(MRInstruments)
  if (nzchar(Sys.getenv("OPENGWAS_JWT")))
    httr::set_config(httr::add_headers(Authorization = paste("Bearer", Sys.getenv("OPENGWAS_JWT"))))
  else log_msg("!! 未设置OPENGWAS_JWT，clump可能失败")
  hs_genes <- unique(c(map_to_human(high_conf$gene), MR_EXTRA_GENES))
  data(gtex_eqtl)
  proxy_tissues <- c("Cells Transformed fibroblasts", "Whole Blood", "Muscle Skeletal")
  exposure_list <- compact(setNames(map(hs_genes, function(g) {
    sub <- subset(gtex_eqtl, gene_name == g & tissue %in% proxy_tissues)
    if (nrow(sub) == 0) return(NULL)
    x <- format_gtex_eqtl(sub); x$F_stat <- (x$beta.exposure / x$se.exposure)^2
    filter(x, F_stat > 10)
  }), hs_genes))
  log_msg("有cis-eQTL工具变量的基因: ", length(exposure_list), "/", length(hs_genes))
  exposure_clumped <- compact(map(exposure_list, function(x)
    tryCatch(clump_data(x, clump_r2 = 0.001, clump_kb = 10000),
             error = function(e) { log_msg("clump失败(", e$message, ")"); NULL })))
  outcome_dat <- extract_outcome_data(snps = unique(unlist(map(exposure_clumped, ~ .x$SNP))),
                                      outcomes = OUTCOME_IDS)
  mr_results <- map(exposure_clumped, function(exp) {
    dat <- harmonise_data(exp, outcome_dat, action = 2)
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0) return(NULL)
    list(mr = mr(dat),
         steiger = tryCatch(directionality_test(dat), error = function(e) NULL),
         harmonised = dat)
  }) %>% compact()
  mr_summary <- map_dfr(mr_results, ~ as.data.frame(.x$mr), .id = "gene")
  write.csv(mr_summary, "output/tables/Table6_MR_results.csv", row.names = FALSE)
  wr <- mr_summary %>% filter(method %in% c("Wald ratio","Inverse variance weighted")) %>%
    select(gene, outcome, method, b, se, pval) %>% arrange(pval)
  print(as.data.frame(wr))
  steiger_df <- map_dfr(mr_results, ~ tryCatch(as.data.frame(.x$steiger), error = function(e) NULL), .id = "gene")
  if (nrow(steiger_df) > 0) steiger_df <- steiger_df %>% distinct(gene, .keep_all = TRUE)
  wr_causal <- wr %>% filter(pval < 0.05) %>% distinct(gene, outcome, .keep_all = TRUE)
  if (nrow(steiger_df) > 0)
    wr_causal <- wr_causal %>%
      left_join(steiger_df %>% select(gene, correct_causal_direction), by = "gene") %>%
      filter(is.na(correct_causal_direction) | correct_causal_direction)
  wr_causal <- wr_causal %>%
    left_join(high_conf %>% mutate(gene_key = toupper(gene)) %>%
                select(gene_key, sc_lfc, slope_inter), by = c("gene" = "gene_key")) %>%
    mutate(direction_match = sign(b) == sign(sc_lfc) |
             (!is.na(slope_inter) & sign(b) == sign(slope_inter)))
  write.csv(wr_causal, "output/tables/Table7_MR_causal_genes.csv", row.names = FALSE)
  log_msg("MR名义阳性(p<0.05，未校正): ", nrow(wr_causal), "个组合")

  wr_all <- mr_summary %>% filter(method == "Wald ratio") %>%
    mutate(lo = b - 1.96 * se, hi = b + 1.96 * se,
           outcome_id = str_extract(outcome, paste(OUTCOME_IDS, collapse = "|")),
           outcome_lab = dplyr::recode(outcome_id, !!!OUTCOME_SHORT, .default = outcome_id))
  if (nrow(wr_all) == 0)
    stop("MR 无 Wald ratio 结果：请检查 clump / extract_outcome_data（需 OPENGWAS_JWT）")
  log_msg("outcome_id → outcome_lab 映射校验: ")
  print(table(wr_all$outcome_id, wr_all$outcome_lab))
  occ_cols <- c("eBMD"="#E64B35","FN-BMD"="#00A087","Osteoporosis"="#7E6148",
                "Fracture"="#4DBBD5","Forearm fx"="#F39B7F","Fracture(self-rpt)"="#8491B4")
  wr_pos <- wr_all %>% filter(pval < 0.05) %>%
    mutate(gene_out = paste(gene, outcome_lab))
  save_fig("Fig9_MR_positives", 10, max(5, 0.30 * nrow(wr_pos) + 2), {
  print(ggplot(wr_pos, aes(b, reorder(gene_out, b), color = outcome_lab)) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) +
        geom_point(size = 2.2) +
        scale_color_manual(values = occ_cols) +
        labs(x = "Wald ratio β (per SD expression, 95% CI)", y = NULL, color = "Outcome",
             title = paste0("cis-MR nominal positives (n = ", nrow(wr_pos), ")")) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))) +
        theme(legend.position = "bottom", legend.box = "horizontal",
              legend.text = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              plot.margin = margin(5, 10, 18, 5),
              plot.title = element_text(size = 11, face = "bold")))
  })
  ## v3.8: position_dodge2 无 height 参数（v3.7 实测报错根源）→ 改 width，并加兜底
  pos_d <- tryCatch(position_dodge2(width = 0.62, preserve = "single"),
                    error = function(e) {
                      log_msg("position_dodge2 失败，回退 position_dodge: ", e$message)
                      position_dodge(width = 0.6) })
  save_fig("FigS10_MR_forest_all", 11, 12, {
  print(ggplot(wr_all, aes(b, gene, color = outcome_lab, group = outcome_lab)) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16,
                       position = pos_d, alpha = 0.85) +
        geom_point(size = 1.8, position = pos_d) +
        scale_color_manual(values = occ_cols) +
        labs(x = "Wald ratio β (per SD expression)", y = NULL, color = "Outcome",
             title = paste0("All ", nrow(wr_all), " cis-MR tests")) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))) +
        theme(legend.position = "bottom", legend.box = "horizontal",
              legend.text = element_text(size = 8),
              axis.text.y = element_text(size = 6.5),
              plot.margin = margin(5, 10, 18, 5),
              plot.title = element_text(size = 11, face = "bold")))
  })
  save(mr_results, mr_summary, wr, wr_causal, wr_all, high_conf, hs_genes,
       file = "output/MR_results.RData")
})

## ---------------- 13. coloc ----------------
library(coloc); library(ieugwasr)
if (file.exists("output/MR_results.RData")) load("output/MR_results.RData")
coloc_genes <- if (exists("wr_causal") && nrow(wr_causal) > 0) {
  dm <- wr_causal$direction_match
  if (!is.null(dm) && any(dm %in% TRUE)) unique(wr_causal$gene[dm %in% TRUE])
  else { log_msg("!! direction_match无TRUE值 → 回退全部阳性基因（请确认MR结果）")
         unique(wr_causal$gene) }
} else character(0)
if (length(coloc_genes) == 0) {
  log_msg("无MR阳性基因 → 跳过coloc（或手动指定 coloc_genes）")
} else {
  COLOC_GWAS <- list(id = OUTCOME_IDS[1], N = 484598, type = "cc", s = 8844/484598)
  EqtlN <- c(Fibroblasts = 483, WholeBlood = 670)
  run_coloc <- function(gene, tissue_name) {
    f <- paste0("coloc_data/eQTL_", gene, "_", tissue_name, ".csv")
    if (!file.exists(f)) { log_msg("缺文件 ", f); return(NULL) }
    eq <- read.csv(f, stringsAsFactors = FALSE) %>%
      mutate(se_e = abs(beta) / pmax(abs(qnorm(pval / 2)), 0.5)) %>%
      filter(is.finite(se_e), se_e > 0)
    gw <- tryCatch(ieugwasr::associations(variants = eq$snp, id = COLOC_GWAS$id, proxies = 0),
                   error = function(e) NULL)
    if (is.null(gw) || nrow(gw) < 10) return(NULL)
    rsid_col <- intersect(c("rsid","variant","snp","name"), names(gw))[1]
    ea_col <- intersect(c("effect_allele","ea"), names(gw))[1]
    oa_col <- intersect(c("other_allele","nea","oa"), names(gw))[1]
    locus <- gw %>%
      transmute(snp = .data[[rsid_col]], beta_g = beta, se_g = se, p_g = p, maf_g = eaf,
                gw_eff = .data[[ea_col]], gw_oth = .data[[oa_col]]) %>%
      inner_join(eq %>% select(snp, beta_e = beta, se_e, eq_ref = ref, eq_alt = alt), by = "snp") %>%
      mutate(flip = case_when(
          toupper(eq_alt) == toupper(gw_eff) & toupper(eq_ref) == toupper(gw_oth) ~ 1L,
          toupper(eq_alt) == toupper(gw_oth) & toupper(eq_ref) == toupper(gw_eff) ~ -1L,
          TRUE ~ 0L)) %>%
      filter(flip != 0L) %>%
      mutate(beta_e = beta_e * flip, MAF = coalesce(maf_g, 0.3))
    if (nrow(locus) < 10) return(NULL)
    d2 <- list(beta = locus$beta_g, varbeta = locus$se_g^2, N = COLOC_GWAS$N,
               type = COLOC_GWAS$type, snp = locus$snp, MAF = locus$MAF)
    if (COLOC_GWAS$type == "cc") d2$s <- COLOC_GWAS$s
    res <- coloc.abf(dataset1 = list(beta = locus$beta_e, varbeta = locus$se_e^2,
                                     N = EqtlN[[tissue_name]], type = "quant",
                                     snp = locus$snp, MAF = locus$MAF),
                     dataset2 = d2)
    sm <- res$summary; names(sm) <- gsub("\\.abf$", "", names(sm))
    tibble(gene = gene, tissue = tissue_name, nsnps = nrow(locus),
           PP.H3 = unname(sm["PP.H3"]), PP.H4 = unname(sm["PP.H4"]),
           colocalized = unname(sm["PP.H4"]) > 0.75)
  }
  coloc_results <- map_dfr(coloc_genes, function(g)
    map_dfr(c("Fibroblasts","WholeBlood"), function(t)
      tryCatch(run_coloc(g, t), error = function(e) NULL)))
  write.csv(coloc_results, "output/tables/TableS11_coloc.csv", row.names = FALSE)
  print(as.data.frame(coloc_results))
  log_msg("coloc完成: 通过PP.H4>0.75的 基因×组织 = ",
          sum(coloc_results$colocalized, na.rm = TRUE))
  save(coloc_results, file = "output/coloc_results.RData")
}

## ---------------- 14. 汇总图 + 存档 ----------------
grp_lbl <- c("Young_Intact"="YI","Young_Fracture"="YF","Aged_Intact"="AI","Aged_Fracture"="AF")
p0a <- mk_donor_panel(sen_intact, ylab = "SenMayo score", title = "(a) Senescent remodeling")
save_fig("Fig0_graphical_abstract", 12, 4.4, {
print(wrap_plots(
  p0a,
  ggplot(pt %>% count(group) %>% mutate(group = factor(group, levels = names(grp_lbl))),
         aes(group, n, fill = group)) + geom_col() + theme(legend.position = "none") +
    scale_x_discrete(labels = grp_lbl, drop = FALSE) +
    scale_fill_manual(values = c("Young_Intact"="#F39B7F","Young_Fracture"="#7CAE00",
                                 "Aged_Intact"="#00BFC4","Aged_Fracture"="#C77CFF")) +
    labs(title = "(b) Abortive activation", x = NULL, y = "Cells",
         caption = "YI/YF/AI/AF = Young/Aged × Intact/Fracture"),
  ggplot(head(sen_comm, 12), aes(reorder(paste(ligand, receptor), prob), prob, fill = age)) +
    geom_col() + coord_flip() + theme(legend.position = "none") +
    scale_fill_manual(values = c("Young"="#00A087","Aged"="#E64B35")) +
    labs(title = "(c) Communication programs", x = NULL, y = "Communication prob."),
  ncol = 3) +
  plot_annotation(title = "Convergent senescence remodeling and abortive regenerative activation in aged periosteum",
                  theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))))
})
save_objs <- c("res_age","deg_age","dds","dds_ix","high_conf","res_ix","gsea_sen",
               "gsea_resp","sen_by_donor","sen_test","auc_val","med_diff","ego_up","ego_dn")
save_objs <- save_objs[vapply(save_objs, exists, logical(1), envir = .GlobalEnv)]
save(list = save_objs, file = "output/key_results_periosteum.RData")
log_msg("key_results_periosteum.RData 已保存")

## v3.8【补丁1/2】为 run_periosteum_strengthening.R 提供输入（带防呆）
if (exists("peri")) {
  saveRDS(peri, "output/peri_seurat.rds")
  log_msg("补丁1: peri_seurat.rds 已保存")
} else log_msg("!! 补丁1 未执行: peri 不存在")
if (exists("res_age_df")) {
  write.csv(res_age_df, "output/tables/TableS0_full_DE_stats.csv", row.names = FALSE)
  log_msg("补丁2: TableS0_full_DE_stats.csv 已保存")
} else log_msg("!! 补丁2 未执行: res_age_df 不存在")

log_msg("========== 流程结束 (periosteum_aging v3.8) ==========")
