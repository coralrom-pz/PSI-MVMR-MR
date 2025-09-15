# ================================
# COLOC Pipeline
# ================================
# Author: Coral Romualdo-Perez
# Date: 2025-09-03
# Notes:
# - Chunks 2-7 can be run together
# - Pseudo layout example used below:
#     project/
#       data/traits/...
#       data/pqtl/UKB_PPP/Hg37/Unformatted/...
#       outputs/coloc/500kb_with_sensitivity/<panel>/...
# ================================


# ================================
# CHUNK 1 — Load packages (silent)
# ================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tibble)
  library(coloc); library(ggplot2); library(ggpubr)
})


# ================================
# CHUNK 2 — Helpers + Standardizers
# ================================
# -- utility: p from beta/se (two-sided)
p_from_beta_se <- function(beta, se) 2 * pnorm(-abs(beta / se))

# -- helper to coerce "chr1"/"1"/"X"/"Y"/"MT" -> integer
.parse_chr_int <- function(x){
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x == "X"]  <- "23"
  x[x == "Y"]  <- "24"
  x[x %in% c("M","MT")] <- "25"   # optional
  suppressWarnings(as.integer(x))
}

# -- header standardizer (robust to many schemas)
standardize_any <- function(df, n_default = NA_real_) {
  nm <- tolower(names(df)); df <- setNames(df, nm); nrows <- nrow(df)
  pick <- function(keys, default = NA_character_) {
    hits <- intersect(keys, nm); if (!length(hits)) return(rep(default, nrows)); as.character(df[[hits[1]]])
  }
  pick_num <- function(keys, default = NA_real_) {
    hits <- intersect(keys, nm); if (!length(hits)) return(rep(default, nrows)); as.numeric(df[[hits[1]]])
  }
  clean_chr <- function(x) { x[x==""] <- NA_character_; x }
  
  # IDs (prefer rsIDs)
  snp <- clean_chr(pick(c("rsid","rsids","hm_rsid")))
  if (all(is.na(snp))) snp <- clean_chr(pick(c("snp","variant_id","hm_variant_id")))
  
  # positions (add chr_hg19/pos_hg19, and parse "chr1")
  chr <- .parse_chr_int(pick(c("chr","chrom","chromosome","hm_chrom","chr_hg19")))
  pos <- suppressWarnings(as.integer(pick(c("pos","position","bp","base_pair_location","genpos",
                                            "hm_pos","pos_hg19"))))
  
  # alleles / effects
  ea  <- toupper(pick(c("effect_allele","ea","a1","allele1","alt","tested_allele",
                        "effect_allele.exposure","effect_allele.outcome","hm_effect_allele")))
  oa  <- toupper(pick(c("other_allele","oa","a2","allele2","ref","ref_allele","non_effect_allele",
                        "other_allele.exposure","other_allele.outcome","hm_other_allele","allele0")))
  beta <- pick_num(c("beta","beta.exposure","beta.outcome","hm_beta","b","effect"))
  se   <- pick_num(c("se","se.exposure","se.outcome","stderr","standard_error","se_beta"))
  
  # frequencies: accept eaf or ImpMAF (UKB/UKB)
  freq_candidates <- nm[grepl("^frq|^freq|^a1freq|^allele1freq|^frq_a_", nm)]
  eaf <- pick_num(c("eaf","maf","a1freq","eaf.exposure","eaf.outcome",
                    "effect_allele_frequency","hm_effect_allele_frequency",
                    "impmaf", freq_candidates))
  
  # N
  nvec <- pick_num(c("n","samplesize","sample_size","samplesize.exposure","samplesize.outcome"))
  nvec[is.na(nvec)] <- n_default
  
  tibble(snp, chr, pos, ea, oa, beta, se, eaf, n = nvec) %>%
    dplyr::filter(!is.na(beta) & !is.na(se) & !is.na(chr) & !is.na(pos) & !is.na(snp)) %>%
    dplyr::mutate(eaf = ifelse(!is.na(eaf) & (eaf < 0 | eaf > 1), NA_real_, eaf))
}


# ================================
# CHUNK 3 — Load & standardize PSI + pQTL tables
# ================================
# -- PSI (loneliness) raw sumstats (pseudo path)
psi_raw <- fread("project/data/traits/Loneliness_EurRel_Imputed.txt")
psi_std <- standardize_any(psi_raw, n_default = 452302) %>%
  select(snp, chr, pos, ea, oa, beta, se, eaf, n)

# -- helper: read → standardize → enforce schema → de-dup per rsID
read_pqtl_std <- function(path, gene, n_default = NA_real_) {
  dt  <- data.table::fread(path, showProgress = FALSE)
  std <- standardize_any(dt, n_default = n_default) %>%
    transmute(
      snp  = as.character(snp),
      chr  = as.integer(chr),
      pos  = as.integer(pos),
      ea   = toupper(as.character(ea)),
      oa   = toupper(as.character(oa)),
      beta = as.numeric(beta),
      se   = as.numeric(se),
      eaf  = as.numeric(eaf),
      n    = as.numeric(n),
      gene = gene
    ) %>%
    arrange(se) %>%                    # if dup rsIDs exist, keep the most precise
    distinct(snp, .keep_all = TRUE)
  std
}

# -- build standardized pQTL in one go (UKB-PPP paths → pseudo)
pqtl_std <- dplyr::bind_rows(
  read_pqtl_std("project/data/pqtl/UKB_PPP/Hg37/Unformatted/IL1RA_with_rsid_allchr_FIXED_hg19_preFormat.txt", "IL1RA"),
  read_pqtl_std("project/data/pqtl/UKB_PPP/Hg37/Unformatted/IL6_with_rsid_allchr_FIXED_hg19_preFormat.txt",   "IL6"),
  read_pqtl_std("project/data/pqtl/UKB_PPP/Hg37/Unformatted/IL6R_with_rsid_allchr_FIXED_hg19_preFormat.txt",  "IL6R"),
  read_pqtl_std("project/data/pqtl/UKB_PPP/Hg37/Unformatted/TNFR1_with_rsid_allchr_FIXED_hg19_preFormat.txt", "TNFR1"),
  read_pqtl_std("project/data/pqtl/UKB_PPP/Hg37/Unformatted/TNFR2_with_rsid_allchr_FIXED_hg19_preFormat.txt", "TNFR2")
) %>%
  group_by(gene) %>%
  mutate(n = ifelse(is.na(n), 35559, n)) %>%  # fallback N per gene (edit to your cohort N)
  ungroup()

# -- gene coordinates (GRCh37)
genes <- tribble(
  ~gene,  ~chr, ~start,      ~end,
  "IL1RA",  2, 113864791, 113891593,
  "IL6",    7,  22765503,  22771621,
  "IL6R",   1, 154377669, 154441926,
  "TNFR1", 12,   6437923,   6451280,
  "TNFR2",  1,  12227060,  12269285
)


# ================================
# CHUNK 4 — Trait-internal alignment helpers
# ================================
# -- align each trait to its minor allele (per tutorial)
align_to_minor <- function(df, prefix){
  if (!("eaf" %in% names(df)) || all(is.na(df$eaf))) {
    df[[paste0("EA_",prefix,"_aligned")]]   <- df$ea
    df[[paste0("NEA_",prefix,"_aligned")]]  <- df$oa
    df[[paste0("MAF_",prefix)]]             <- NA_real_
    df[[paste0("BETA_",prefix,"_aligned")]] <- df$beta
    return(df)
  }
  df[[paste0("EA_",prefix,"_aligned")]]   <- ifelse(df$eaf < 0.5, df$ea, df$oa)
  df[[paste0("NEA_",prefix,"_aligned")]]  <- ifelse(df$eaf < 0.5, df$oa, df$ea)
  df[[paste0("MAF_",prefix)]]             <- ifelse(df$eaf < 0.5, df$eaf, 1 - df$eaf)
  df[[paste0("BETA_",prefix,"_aligned")]] <- ifelse(df$eaf < 0.5, df$beta, -df$beta)
  df
}


# ================================
# CHUNK 5 — COLOC runner (per gene)
# ================================
run_coloc_gene_supervisor_style <- function(gene_row, window_kb = 500,
                                            out_dir = "project/outputs/coloc/500kb_with_sensitivity/UKB",
                                            psi_default_N = 452302,
                                            pqtl_default_N = 21758){
  g <- gene_row$gene; chr <- gene_row$chr
  w_start <- gene_row$start - window_kb*1000
  w_end   <- gene_row$end   + window_kb*1000
  message(sprintf("[COLOC] %s chr%s:%d-%d (±%dkb)", g, chr, w_start, w_end, window_kb))
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "per_gene_tables"), showWarnings = FALSE)
  dir.create(file.path(out_dir, "plots"),           showWarnings = FALSE)
  
  # A) Cytokine cis window (pQTL)
  pqtl_cis <- pqtl_std %>% filter(gene==g, chr==!!chr, pos>=w_start, pos<=w_end) %>%
    arrange(se) %>% distinct(snp, .keep_all = TRUE)
  # B) PSI by chromosome only
  psi_chr  <- psi_std  %>% filter(chr==!!chr) %>%
    arrange(se) %>% distinct(snp, .keep_all = TRUE)
  
  if(nrow(pqtl_cis) < 50 || nrow(psi_chr) < 50)
    return(tibble(gene=g, chr=chr, n_common=0, note="too_few_snps_prejoin",
                  PP.H4=NA_real_, PP.H3=NA_real_))
  
  # Align within trait, compute p
  pqtl_cis <- align_to_minor(pqtl_cis, "pqtl") %>%
    mutate(P_pqtl = p_from_beta_se(beta, se)) %>%
    select(snp, chr, pos_pqtl=pos,
           beta_pqtl_aligned=BETA_pqtl_aligned, se_pqtl=se,
           maf_pqtl=MAF_pqtl, n_pqtl=n, P_pqtl)
  
  psi_chr  <- align_to_minor(psi_chr,  "psi") %>%
    mutate(P_psi = p_from_beta_se(beta, se)) %>%
    select(snp, chr, pos_psi=pos,
           beta_psi_aligned=BETA_psi_aligned, se_psi=se,
           maf_psi=MAF_psi, n_psi=n, P_psi)
  
  # Merge by rsID + chr
  d <- inner_join(pqtl_cis, psi_chr, by=c("snp","chr"))
  
  # Keep only SNPs still inside cis window by pQTL position (safety)
  d <- d %>% filter(pos_pqtl >= w_start, pos_pqtl <= w_end)
  
  if (nrow(d) < 30)
    return(tibble(gene=g, chr=chr, n_common=nrow(d), note="too_few_after_merge",
                  PP.H4=NA_real_, PP.H3=NA_real_))
  
  # coloc.abf (quant vs quant)
  N1 <- ifelse(all(is.na(d$n_psi)),  psi_default_N, stats::median(d$n_psi,  na.rm=TRUE))
  N2 <- ifelse(all(is.na(d$n_pqtl)), pqtl_default_N, stats::median(d$n_pqtl, na.rm=TRUE))
  
  D1 <- list(type="quant",
             beta = d$beta_psi_aligned,
             varbeta = d$se_psi^2,
             N = N1,
             snp = d$snp,
             sdY = 1)
  if (!anyNA(d$maf_psi)) D1$MAF <- d$maf_psi
  
  D2 <- list(type="quant",
             beta = d$beta_pqtl_aligned,
             varbeta = d$se_pqtl^2,
             N = N2,
             snp = d$snp,
             sdY = 1)
  if (!anyNA(d$maf_pqtl)) D2$MAF <- d$maf_pqtl
  
  coloc_obj <- coloc.abf(D1, D2, p1=1e-4, p2=1e-4, p12=1e-5)
  
  # Save: coloc object, merged table, plots
  rds_path  <- file.path(out_dir, sprintf("%s_chr%s_%d_%d_coloc.rds", g, chr, w_start, w_end))
  saveRDS(list(coloc = coloc_obj, D1 = D1, D2 = D2), rds_path)
  
  out_tab <- d %>% transmute(snp, chr,
                             pos = pos_pqtl,
                             beta_psi = beta_psi_aligned, se_psi, P_psi,
                             beta_pqtl = beta_pqtl_aligned, se_pqtl, P_pqtl,
                             maf_psi, maf_pqtl)
  
  data_path <- file.path(out_dir, "per_gene_tables",
                         sprintf("%s_chr%s_%d_%d_merged.tsv", g, chr, w_start, w_end))
  data.table::fwrite(out_tab, data_path, sep="\t")
  
  p1 <- ggplot(out_tab, aes(x = pos, y = -log10(P_psi))) +
    geom_point() + theme_bw() +
    labs(title = "PSI", x = "", y = expression(-log[10](p))) +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 <- ggplot(out_tab, aes(x = pos, y = -log10(P_pqtl))) +
    geom_point() + theme_bw() +
    labs(title = paste0(g, " pQTL"), x = paste0("Chromosome ", chr, " Position"),
         y = expression(-log[10](p))) +
    theme(plot.title = element_text(hjust = 0.5))
  
  fig <- ggpubr::ggarrange(p1, p2, heights = c(1,1), nrow = 2, ncol = 1, align = "hv")
  png_path <- file.path(out_dir, "plots",
                        sprintf("%s_chr%s_%d_%d.png", g, chr, w_start, w_end))
  ggplot2::ggsave(png_path, plot = fig, width = 8, height = 6, dpi = 300)
  
  tibble(gene=g, chr=chr, start=w_start, end=w_end,
         n_common = nrow(d),
         PP.H0 = coloc_obj$summary["PP.H0.abf"],
         PP.H1 = coloc_obj$summary["PP.H1.abf"],
         PP.H2 = coloc_obj$summary["PP.H2.abf"],
         PP.H3 = coloc_obj$summary["PP.H3.abf"],
         PP.H4 = coloc_obj$summary["PP.H4.abf"],
         coloc_rds = rds_path,
         table_tsv = data_path,
         plot_png  = png_path)
}


# ================================
# CHUNK 6 — Run across cytokines + save per-panel summary
# ================================
coloc_out_dir <- "project/outputs/coloc/500kb_with_sensitivity/UKB"

res <- genes %>%
  rowwise() %>%
  do(run_coloc_gene_supervisor_style(., window_kb = 500, out_dir = coloc_out_dir)) %>%
  ungroup() %>% as.data.frame()

data.table::fwrite(res, file.path(coloc_out_dir, "coloc_summary.tsv"), sep = "\t")
print(res[order(-res$PP.H4), c("gene","PP.H4","n_common","plot_png")])


# ================================
# CHUNK 7 — Sensitivity analysis (p12 sweep, plots)
# ================================
base_dir <- coloc_out_dir
sens_dir <- file.path(base_dir, "sensitivity")
dir.create(sens_dir, showWarnings = FALSE, recursive = TRUE)

# choose any rule(s) you want reported
rules <- c("H4 > 0.8")  # add "H4 > 0.5", "H4 > 3*H3 & H0 < 0.1", etc.

sens_rows <- lapply(seq_len(nrow(res)), function(i){
  gene <- res$gene[i]
  rpath <- res$coloc_rds[i]
  obj   <- readRDS(rpath)
  
  # support both save formats
  coloc_obj <- if (!is.null(obj$coloc)) obj$coloc else obj
  D1        <- if (!is.null(obj$D1)) obj$D1 else NULL
  D2        <- if (!is.null(obj$D2)) obj$D2 else NULL
  
  # 1) built-in sensitivity() plot for the first rule (for convenience)
  png(file.path(sens_dir, sprintf("%s_sensitivity_H4gt0.8.png", gene)),
      width = 900, height = 600)
  print(coloc::sensitivity(coloc_obj, rule = rules[1]))
  dev.off()
  
  # 2) sweep a p12 grid (only if D1/D2 were saved)
  p12_grid <- 10^seq(-8, -2, by = 0.25)
  grid_df  <- NULL
  min_p12_H4gt05 <- NA_real_
  min_p12_H4gt08 <- NA_real_
  max_H4         <- as.numeric(coloc_obj$summary["PP.H4.abf"])
  
  if (!is.null(D1) && !is.null(D2)) {
    grid_df <- data.frame(p12 = p12_grid,
                          PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_,
                          PP.H3 = NA_real_, PP.H4 = NA_real_)
    for (j in seq_along(p12_grid)) {
      rr <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = p12_grid[j])
      post <- as.numeric(rr$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")])
      grid_df[j, 2:6] <- post
    }
    data.table::fwrite(grid_df,
                       file.path(sens_dir, sprintf("%s_p12_grid.tsv", gene)),
                       sep = "\t")
    
    if (any(grid_df$PP.H4 > 0.5)) min_p12_H4gt05 <- min(grid_df$p12[grid_df$PP.H4 > 0.5])
    if (any(grid_df$PP.H4 > 0.8)) min_p12_H4gt08 <- min(grid_df$p12[grid_df$PP.H4 > 0.8])
    max_H4 <- max(grid_df$PP.H4, na.rm = TRUE)
  }
  
  data.frame(
    gene = gene,
    H4_at_default = as.numeric(coloc_obj$summary["PP.H4.abf"]),
    max_H4_across_grid = max_H4,
    min_p12_for_H4_gt_0_5 = min_p12_H4gt05,
    min_p12_for_H4_gt_0_8 = min_p12_H4gt08,
    sens_plot = file.path(sens_dir, sprintf("%s_sensitivity_H4gt0.8.png", gene))
  )
})

sens_summary <- do.call(rbind, sens_rows)
data.table::fwrite(sens_summary, file.path(sens_dir, "sensitivity_summary.tsv"), sep = "\t")
sens_summary


# ================================
# CHUNK 8 — Results combiner across panels
# ================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(tidyr); library(purrr)
})

# --- point these at your output folders (pseudo paths) ---
panel_dirs <- tribble(
  ~panel,        ~dir,
  "Folkersen",   "project/outputs/coloc/500kb_with_sensitivity/Folkersen",
  "Ferkingstad", "project/outputs/coloc/500kb_with_sensitivity/Ferkingstad",
  "UKB_PPP",     "project/outputs/coloc/500kb_with_sensitivity/UKB"
)

read_panel <- function(panel, base){
  summ_p <- file.path(base, "coloc_summary.tsv")
  sens_p <- file.path(base, "sensitivity", "sensitivity_summary.tsv")
  
  if (!file.exists(summ_p)) stop("Missing coloc_summary.tsv in: ", base)
  cs  <- fread(summ_p, showProgress = FALSE)
  cs$panel <- panel
  
  ss  <- if (file.exists(sens_p)) fread(sens_p, showProgress = FALSE) else
    data.table(gene=character(), H4_at_default=numeric(), max_H4_across_grid=numeric(),
               min_p12_for_H4_gt_0_5=numeric(), min_p12_for_H4_gt_0_8=numeric(), sens_plot=character())
  out <- cs %>%
    left_join(ss, by = "gene") %>%
    relocate(panel, gene)
  
  # add quick classification flags
  out <- out %>%
    mutate(
      PP.H0 = as.numeric(PP.H0), PP.H1 = as.numeric(PP.H1),
      PP.H2 = as.numeric(PP.H2), PP.H3 = as.numeric(PP.H3),
      PP.H4 = as.numeric(PP.H4),
      max_H4_across_grid   = as.numeric(max_H4_across_grid),
      min_p12_for_H4_gt_0_5= as.numeric(min_p12_for_H4_gt_0_5),
      min_p12_for_H4_gt_0_8= as.numeric(min_p12_for_H4_gt_0_8),
      top_h = pmax(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, na.rm = TRUE),
      call = dplyr::case_when(
        PP.H4 >= 0.80                                   ~ "Robust colocalization (H4)",
        PP.H3 >= 0.50                                   ~ "Likely LD confounding (H3)",
        PP.H4 >= 0.50 & PP.H4 < 0.80                    ~ "Suggestive coloc (sensitive)",
        PP.H2 >= 0.50 & PP.H3 < 0.20 & PP.H4 < 0.20     ~ "Exposure-only (pQTL-only)",
        PP.H1 >= 0.50 & PP.H3 < 0.20 & PP.H4 < 0.20     ~ "Outcome-only",
        PP.H0 >= 0.50                                   ~ "No regional signal",
        TRUE                                            ~ "Mixed/uncertain"
      ),
      robustness = dplyr::case_when(
        PP.H4 >= 0.80                                                           ~ "Robust @ defaults",
        PP.H4 < 0.80 & !is.na(min_p12_for_H4_gt_0_8) & min_p12_for_H4_gt_0_8 <= 1e-4
        ~ "Becomes robust if p12≤1e-4",
        PP.H4 < 0.50 & !is.na(max_H4_across_grid) & max_H4_across_grid < 0.50   ~ "No coloc across p12 grid",
        TRUE                                                                    ~ "Check sensitivity plot"
      )
    )
  out
}

all_panels <- panel_dirs %>% pmap_dfr(~ read_panel(..1, ..2))

# Optional: compute z-score concordance if merged tables exist
all_panels <- all_panels %>%
  rowwise() %>%
  mutate(z_corr = {
    dir <- panel_dirs$dir[match(panel, panel_dirs$panel)]
    files <- list.files(file.path(dir, "per_gene_tables"),
                        pattern = paste0("^", gene, "_.*_merged\\.tsv$"),
                        full.names = TRUE)
    if (!length(files)) NA_real_ else {
      tab <- data.table::fread(files[1], showProgress = FALSE)
      if (!all(c("beta_psi","se_psi","beta_pqtl","se_pqtl") %in% names(tab))) NA_real_ else {
        cor(tab$beta_psi/tab$se_psi, tab$beta_pqtl/tab$se_pqtl, use="complete.obs")
      }
    }
  }) %>%
  ungroup()

# Save a tidy comparison table
dir.create("project/outputs/coloc/summary_across_panels", showWarnings = FALSE, recursive = TRUE)
fwrite(all_panels, "project/outputs/coloc/summary_across_panels/coloc_across_panels.tsv", sep = "\t")

# Show a concise view
all_panels %>%
  arrange(gene, panel) %>%
  select(panel, gene, n_common, PP.H2, PP.H3, PP.H4, max_H4_across_grid, robustness, call)