#!/usr/bin/env Rscript
## ============================================================
## run_coloc_headline4_eBMD.R  (v1.0)
## ------------------------------------------------------------
## 目的：为 4 个 Bonferroni 显著基因（DOCK9 / SPP1 / TMTC2 / APLNR）
##       补做 cis-eQTL ↔ 结局 GWAS 共定位，填补 Table S26/S27 空白。
##
## 相对原流程的改动（对应审稿隐患）：
##   ① 原 §13 只跑 direction_match==TRUE 的基因 → 本脚本直接指定 4 个头部分子；
##   ② 原共定位结局是 Fractures(GCST90038703) → 本脚本主分析改为
##      heel eBMD(ebi-a-GCST90029004, quant, N=583314)，与 MR 的 Bonferroni
##      结论同结局；Fractures 保留为敏感性分析；
##   ③ coloc.susie 沿用 strengthening §11 策略：尝试 ieugwasr::ld_matrix
##      获取 EUR LD，失败记 NA（与正文 "in-sample LD unavailable" 一致）。
##
## 复用与缓存：
##   - GTEx v8 signif_pairs 与 variant lookup 复用 coloc_data/_cache/，
##     已下载过则秒过；缺失才自动下载（文件很大，需耐心/代理）。
##   - 已存在的 coloc_data/eQTL_<GENE>_<组织>.csv 默认不重提
##     （REFRESH_EQTL <- TRUE 可强制重提）。
##
## 输出（追加格式与 Table S26 / S27 完全一致）：
##   coloc_data/eQTL_<GENE>_<Fibroblasts|WholeBlood>.csv
##   output/tables/coloc_headline4_eBMD.csv          -> 追加进 TableS26_coloc.csv
##   output/tables/coloc_headline4_eBMD_susie.csv    -> 追加进 TableS27_coloc_susie.csv
##   output/tables/coloc_headline4_fractures.csv     （敏感性，可选）
##   output/tables/coloc_headline4_fractures_susie.csv（敏感性，可选）
##
## 运行前提：工作目录 = 主分析项目根（与 run_periosteum_aging20.R 相同，
##            即包含 output/ 和 coloc_data/ 的目录）
## 运行方式：Rscript run_coloc_headline4_eBMD.R
## 内存提示：读取 variant lookup 约需 8–16 GB 空闲内存（与原提取脚本相同）
## ============================================================

## ---------------- 0. 配置 ----------------
GENES            <- c("DOCK9", "SPP1", "TMTC2", "APLNR")
TISSUES          <- c("Fibroblasts", "WholeBlood")
RUN_SENSITIVITY  <- TRUE            # 是否另跑 Fractures 敏感性分析
REFRESH_EQTL     <- FALSE           # TRUE = 无视已存在的 eQTL csv 强制重提
MIN_SNP_ABF      <- 10              # 与主脚本 §13 一致
MIN_SNP_SUSIE    <- 50              # 与 strengthening §11 一致
set.seed(123)

OUT_EBMD <- list(id = "ebi-a-GCST90029004", N = 583314, type = "quant", s = NULL)
OUT_FRAC <- list(id = "ebi-a-GCST90038703", N = 484598, type = "cc",    s = 8844/484598)
EqtlN <- c(Fibroblasts = 483, WholeBlood = 670)   # GTEx v8 组织样本量

suppressPackageStartupMessages({
  library(org.Hs.eg.db)   # 先加载：它的 S4 select 会被下面的 tidyverse 遮蔽（顺序不可换）
  library(data.table); library(tidyverse)
  library(coloc); library(ieugwasr)
})
dir.create("coloc_data/_cache", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables",    showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

## ---------------- 1. 基因 → Ensembl 映射 ----------------
ens <- suppressMessages(mapIds(org.Hs.eg.db, keys = GENES, column = "ENSEMBL",
                               keytype = "SYMBOL", multiVals = "first"))
missing_sym <- GENES[is.na(ens)]
if (length(missing_sym) > 0) stop("Ensembl 映射失败: ", paste(missing_sym, collapse = ", "))
ens <- ens[!is.na(ens)]
sym_map <- setNames(names(ens), unname(ens))
log_msg("基因清单(", length(GENES), "): ", paste(GENES, collapse = ", "))
log_msg("Ensembl: ", paste(paste(names(ens), unname(ens), sep = "="), collapse = ", "))

## ---------------- 2. GTEx v8 缓存解析（复用 fetch_coloc_eqtl.R 逻辑） ----------------
URLS <- c(
  Fibroblasts = "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL/Cells_Cultured_fibroblasts.v8.signif_variant_gene_pairs.txt.gz",
  WholeBlood  = "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"
)
LOOKUP_URL <- "https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"

resolve_cache <- function(url, alias) {
  p_orig <- file.path("coloc_data/_cache", basename(url))
  p_al   <- file.path("coloc_data/_cache", alias)
  for (p in c(p_al, p_orig)) if (file.exists(p) && file.size(p) > 1e6) {
    log_msg("缓存命中: ", basename(p)); return(p)
  }
  log_msg("下载: ", basename(url), " （文件很大，可能需要代理）...")
  code <- tryCatch(download.file(url, p_al, mode = "wb", method = "libcurl", timeout = 7200),
                   error = function(e) { message("  错误: ", e$message); 1 })
  if (code != 0) stop("下载失败: ", url,
    "\n  → 用浏览器+代理手动下载，放入 coloc_data/_cache/（保持.gz，文件名 ",
    alias, " 或原名）后重跑")
  p_al
}

need_extract <- REFRESH_EQTL || any(!file.exists(sprintf(
  "coloc_data/eQTL_%s_%s.csv", rep(GENES, each = length(TISSUES)), TISSUES)))

if (need_extract) {
  lk_path <- resolve_cache(LOOKUP_URL, "gtex_v8_variant_lookup.txt.gz")
  lk_hd <- names(fread(lk_path, nrows = 0))
  vid_col <- grep("variant_id", lk_hd, value = TRUE)[1]
  rs_col  <- grep("rs_id|rsid", lk_hd, value = TRUE)[1]
  ref_col <- lk_hd[tolower(lk_hd) == "ref"][1]
  alt_col <- lk_hd[tolower(lk_hd) == "alt"][1]
  sel <- unique(c(vid_col, rs_col, ref_col, alt_col)); sel <- sel[!is.na(sel)]
  log_msg("读取 variant lookup（", paste(sel, collapse = ", "), "，约 2–4 分钟，需 8–16 GB 内存）...")
  lk_dt <- fread(lk_path, select = sel)
  setnames(lk_dt, vid_col, "variant_id")
  if (!is.na(rs_col) && rs_col != "variant_id") setnames(lk_dt, rs_col, "rs_id")
  log_msg("lookup 行数: ", nrow(lk_dt))

  for (tissue in TISSUES) {
    p <- resolve_cache(URLS[tissue], paste0(tissue, ".v8.signif_pairs.txt.gz"))
    hd <- names(fread(p, nrows = 0))
    need <- intersect(c("variant_id", "gene_id", "pval_nominal", "slope"), hd)
    if (length(need) < 4) stop(tissue, " 列名不符预期: ", paste(hd, collapse = ", "))
    log_msg("读取 ", tissue, " signif_pairs ...")
    dt <- fread(p, select = need)
    dt[, gene_base := sub("\\..*$", "", gene_id)]
    dt <- dt[gene_base %in% unname(ens)]
    log_msg("  目标基因行数: ", nrow(dt))
    if (nrow(dt) == 0) next

    m <- merge(dt, lk_dt, by = "variant_id", all.x = TRUE)
    if (!"rs_id" %in% names(m) || all(is.na(m$rs_id))) {
      parts <- tstrsplit(m$variant_id, "_")
      m[, rs_id := paste0("chr", parts[[2]], ":", parts[[3]])]  ## 兜底 chr:pos（OpenGWAS 可解析）
    }
    if (!"ref" %in% names(m) || all(is.na(m$ref))) {
      parts <- tstrsplit(m$variant_id, "_")
      m[, `:=`(ref = parts[[4]], alt = parts[[5]])]
    }
    m <- m[rs_id != "." & !is.na(rs_id)]

    for (g in unique(m$gene_base)) {
      gsym <- sym_map[[g]]
      sub2 <- m[gene_base == g, .(snp = rs_id, beta = slope, pval = pval_nominal, ref, alt)]
      sub2 <- sub2[complete.cases(sub2) & is.finite(beta) & pval > 0 & ref != alt]
      if (nrow(sub2) < MIN_SNP_ABF) {
        log_msg("  !! ", gsym, " (", tissue, ") 仅 ", nrow(sub2), " SNP（<", MIN_SNP_ABF, "）→ 写出但预计无法共定位")
      }
      out <- sprintf("coloc_data/eQTL_%s_%s.csv", gsym, tissue)
      fwrite(unique(sub2), out)
      log_msg("  写出: ", out, " (", nrow(unique(sub2)), " SNP)")
    }
  }
  rm(lk_dt, dt, m); gc(reset = TRUE)
} else {
  log_msg("4 个基因的 eQTL 文件均已存在 → 跳过提取（REFRESH_EQTL=TRUE 可强制重提）")
}

## ---------------- 3. 共定位（abf 必跑 + susie 尝试） ----------------
run_coloc_one <- function(gene, tissue, oc) {
  f <- sprintf("coloc_data/eQTL_%s_%s.csv", gene, tissue)
  if (!file.exists(f)) { log_msg("  缺 eQTL 文件: ", f, " → 跳过"); return(list(abf = NULL, susie = NULL)) }

  eq <- read.csv(f, stringsAsFactors = FALSE) %>%
    mutate(se_e = abs(beta) / pmax(abs(qnorm(pval / 2)), 0.5)) %>%
    filter(is.finite(se_e), se_e > 0)

  gw <- tryCatch(ieugwasr::associations(variants = eq$snp, id = oc$id, proxies = 0),
                 error = function(e) { log_msg("  GWAS 查询失败(", oc$id, "): ", e$message); NULL })
  if (is.null(gw) || nrow(gw) < MIN_SNP_ABF) {
    log_msg("  ", gene, "/", tissue, ": GWAS 区域 SNP 不足 (", if (is.null(gw)) 0 else nrow(gw), ") → 跳过")
    return(list(abf = NULL, susie = NULL))
  }
  rsid_col <- intersect(c("rsid", "variant", "snp", "name"), names(gw))[1]
  ea_col   <- intersect(c("effect_allele", "ea"), names(gw))[1]
  oa_col   <- intersect(c("other_allele", "nea", "oa"), names(gw))[1]
  locus <- gw %>%
    transmute(snp = .data[[rsid_col]], beta_g = beta, se_g = se,
              maf_g = eaf, gw_eff = .data[[ea_col]], gw_oth = .data[[oa_col]]) %>%
    inner_join(eq %>% dplyr::select(snp, beta_e = beta, se_e, eq_ref = ref, eq_alt = alt), by = "snp") %>%
    mutate(flip = case_when(
      toupper(eq_alt) == toupper(gw_eff) & toupper(eq_ref) == toupper(gw_oth) ~ 1L,
      toupper(eq_alt) == toupper(gw_oth) & toupper(eq_ref) == toupper(gw_eff) ~ -1L,
      TRUE ~ 0L)) %>%
    filter(flip != 0L) %>%
    mutate(beta_e = beta_e * flip, MAF = coalesce(maf_g, 0.3))
  if (nrow(locus) < MIN_SNP_ABF) {
    log_msg("  ", gene, "/", tissue, ": 回交后 SNP 不足 (", nrow(locus), ") → 跳过")
    return(list(abf = NULL, susie = NULL))
  }

  ## coloc.abf（单信号，主结果）
  d2 <- list(beta = locus$beta_g, varbeta = locus$se_g^2, N = oc$N,
             type = oc$type, snp = locus$snp, MAF = locus$MAF)
  if (oc$type == "cc") d2$s <- oc$s
  abf_res <- tryCatch(coloc::coloc.abf(
    dataset1 = list(beta = locus$beta_e, varbeta = locus$se_e^2, N = EqtlN[[tissue]],
                    type = "quant", snp = locus$snp, MAF = locus$MAF),
    dataset2 = d2)$summary, error = function(e) { log_msg("  coloc.abf 失败: ", e$message); NULL })

  abf_row <- NULL
  if (!is.null(abf_res)) {
    names(abf_res) <- gsub("\\.abf$", "", names(abf_res))
    abf_row <- tibble(gene = gene, tissue = tissue, nsnps = nrow(locus),
                      PP.H3 = unname(abf_res["PP.H3"]), PP.H4 = unname(abf_res["PP.H4"]),
                      colocalized = unname(abf_res["PP.H4"]) > 0.75)
  }

  ## coloc.susie（多信号；LD 获取失败则全 NA，与主文表述一致）
  sus_row <- tibble(gene = gene, tissue = tissue, nsnps = nrow(locus),
                    PP_H4_abf = if (!is.null(abf_res)) unname(abf_res["PP.H4"]) else NA_real_,
                    PP_H4_susie = NA_real_, n_signals_eqtl = NA_integer_,
                    n_signals_gwas = NA_integer_, susie_pass = FALSE)
  LD_mat <- NULL
  if (nrow(locus) >= MIN_SNP_SUSIE) {
    LD_mat <- tryCatch(ieugwasr::ld_matrix(locus$snp, pop = "EUR"),
                       error = function(e) { log_msg("  LD 获取失败(", gene, ",", tissue, "): ", e$message); NULL })
    if (!is.null(LD_mat)) {
      cs <- intersect(locus$snp, rownames(LD_mat))
      if (length(cs) >= MIN_SNP_SUSIE) {
        LD_mat <- LD_mat[cs, cs, drop = FALSE]
        locus  <- locus[match(cs, locus$snp), , drop = FALSE]
      } else LD_mat <- NULL
    }
  }
  if (!is.null(LD_mat)) {
    D1 <- list(beta = locus$beta_e, varbeta = locus$se_e^2, N = EqtlN[[tissue]],
               type = "quant", snp = locus$snp, MAF = locus$MAF, LD = LD_mat)
    D2 <- list(beta = locus$beta_g, varbeta = locus$se_g^2, N = oc$N, type = oc$type,
               snp = locus$snp, MAF = locus$MAF, LD = LD_mat)
    if (oc$type == "cc") D2$s <- oc$s
    S1 <- tryCatch(coloc::runsusie(D1), error = function(e) { log_msg("  runsusie eQTL 失败: ", e$message); NULL })
    S2 <- tryCatch(coloc::runsusie(D2), error = function(e) { log_msg("  runsusie GWAS 失败: ", e$message); NULL })
    if (!is.null(S1) && !is.null(S2) && length(S1@pip) > 0 && length(S2@pip) > 0) {
      sus_res <- tryCatch(coloc::coloc.susie(S1, S2)$summary, error = function(e) NULL)
      if (!is.null(sus_res)) {
        sus_row$PP_H4_susie   <- unname(sus_res["PP.H4.abf"])
        sus_row$n_signals_eqtl  <- length(S1@pip[S1@pip > 0.5])
        sus_row$n_signals_gwas  <- length(S2@pip[S2@pip > 0.5])
        sus_row$susie_pass      <- unname(sus_res["PP.H4.abf"]) > 0.75
      }
    }
  } else {
    log_msg("  ", gene, "/", tissue, ": susie 跳过（LD 不可用或无足够 SNP）→ 记 NA")
  }
  list(abf = abf_row, susie = sus_row)
}

run_outcome <- function(oc, tag) {
  log_msg("==== 共定位：结局 ", oc$id, " (", tag, ") ====")
  res <- purrr::map(GENES, function(g)
    purrr::map(TISSUES, function(t) run_coloc_one(g, t, oc))) %>% unlist(recursive = FALSE)
  abf  <- bind_rows(lapply(res, `[[`, "abf"))
  sus  <- bind_rows(lapply(res, `[[`, "susie"))
  write.csv(abf, sprintf("output/tables/coloc_headline4_%s.csv", tag),        row.names = FALSE)
  write.csv(sus, sprintf("output/tables/coloc_headline4_%s_susie.csv", tag), row.names = FALSE)
  log_msg(tag, " 完成。coloc.abf 结果："); print(as.data.frame(abf))
  if (nrow(abf) > 0 && any(abf$colocalized, na.rm = TRUE))
    warning("!! 出现 PP.H4 > 0.75 的基因-组织对 → 正文与 Discussion 措辞需重写，请立即告知")
  list(abf = abf, susie = sus)
}

main_ebmd <- run_outcome(OUT_EBMD, "eBMD")
if (RUN_SENSITIVITY) run_outcome(OUT_FRAC, "fractures")

## ---------------- 4. 收尾提示 ----------------
log_msg("全部完成。下一步：")
log_msg("  1) 检查 output/tables/coloc_headline4_eBMD.csv 与 _susie.csv；")
log_msg("  2) 把行追加进 TableS26_coloc.csv / TableS27_coloc_susie.csv（保持列一致），")
log_msg("     并同步更新补充稿 Table S26/S27 标题中的对数；")
log_msg("  3) 把两个 CSV 发回，更新正文 P38 定稿措辞并撤掉红色批注。")
sessionInfo()
