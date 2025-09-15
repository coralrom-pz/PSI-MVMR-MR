# ================================
# UKB-PPP pQTL: Liftover + Format (hg38 -> hg19) + cis windows + MR exposure
# ================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(tidyr)
  library(GenomicRanges); library(rtracklayer); library(TwoSampleMR); library(ieugwasr)
  library(genetics.binaRies)  # for PLINK path
})

# ---- paths (edit) ----
IN_DIR     <- "/Volumes/Marohu/UKB-PPP/Annotated Cytokines_Combined/Fixed/CleanedInputs"  # hg38 input files ( *_CLEANED.txt )
OUT_DIR    <- "/Volumes/Marohu/UKB-PPP/Annotated Cytokines_Combined/Fixed/Hg37"           # all hg19 outputs
CHAIN_PATH <- "/Volumes/Marohu/hg38ToHg19.over.chain.gz"                                   # must be .gz
REF_BFILE  <- "/Users/coralrom/Desktop/LD Clumping Data/EUR"                               # 1KG EUR (hg19) prefix
PLINK_BIN  <- genetics.binaRies::get_plink_binary()

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# CHUNK 1 — Helpers
# =========================
to_num  <- function(x) suppressWarnings(as.numeric(x))
is_rsid <- function(x) !is.na(x) & grepl("^rs[0-9]+$", x)

# Standardize a UKB-PPP header (hg38); compute p from LOG10P; split multi-RSIDs
standardize_ukb_v38 <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, fill = TRUE, data.table = FALSE, check.names = FALSE)
  names(dt) <- gsub("\\.", "", names(dt))
  
  need <- c("CHROM","GENPOS","RSID","ALLELE1","ALLELE0","BETA","SE","LOG10P","N")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing columns in ", basename(path), ": ", paste(miss, collapse = ", "))
  
  # compute p from LOG10P; clamp to [1e-300, 1]
  Pval <- 10^(-to_num(dt$LOG10P))
  Pval <- pmin(pmax(Pval, 1e-300), 1)
  
  # split multi-rsids like "rs1;rs2"
  dt2 <- tibble(
    chrom = dt$CHROM,
    pos   = dt$GENPOS,
    rsids = ifelse(dt$RSID == "", NA, dt$RSID),
    ea    = toupper(dt$ALLELE1),
    oa    = toupper(dt$ALLELE0),
    beta  = to_num(dt$BETA),
    se    = to_num(dt$SE),
    pval  = Pval,
    n     = to_num(dt$N)
  ) %>%
    tidyr::separate_rows(rsids, sep = ";")
  
  tibble(
    chrom38 = paste0("chr", as.integer(dt2$chrom)),  # GRanges expects "chrN"
    pos38   = as.integer(dt2$pos),
    rsids   = dt2$rsids,
    ea      = dt2$ea,
    oa      = dt2$oa,
    beta    = dt2$beta,
    se      = dt2$se,
    pval    = dt2$pval,
    n       = dt2$n,
    eaf     = NA_real_        # unknown in these inputs; OK for MR/coloc
  )
}

# Liftover hg38 -> hg19; keep only 1:1 mappings
liftover_to_hg19 <- function(df38, chain_path = CHAIN_PATH) {
  gr38 <- GRanges(seqnames = df38$chrom38, ranges = IRanges(df38$pos38, df38$pos38))
  ch   <- import.chain(chain_path)
  lo   <- liftOver(gr38, ch)
  mapn <- lengths(lo)
  
  keep <- which(mapn == 1L)
  out  <- df38[keep, , drop = FALSE]
  lo1  <- unlist(lo[keep])
  
  out$chr_hg19 <- as.character(seqnames(lo1))
  out$pos_hg19 <- start(lo1)
  
  list(
    mapped  = out,
    summary = list(total = nrow(df38), onetoone = length(keep),
                   multi = sum(mapn > 1L), unmapped = sum(mapn == 0L))
  )
}

# Build “preFormat” frame (prefer rsID; else synthesize chr:pos:EA:OA)
build_preformat <- function(df_mapped) {
  snp <- ifelse(is_rsid(df_mapped$rsids) & df_mapped$rsids != "NA",
                df_mapped$rsids,
                paste0(gsub("^chr","", df_mapped$chr_hg19), ":", df_mapped$pos_hg19, ":",
                       df_mapped$ea, ":", df_mapped$oa))
  tibble(
    SNP           = snp,
    effect_allele = df_mapped$ea,
    other_allele  = df_mapped$oa,
    beta          = df_mapped$beta,
    se            = df_mapped$se,
    pval          = pmin(pmax(df_mapped$pval, 1e-300), 1),
    eaf           = df_mapped$eaf,
    samplesize    = df_mapped$n,
    chr_hg19      = gsub("^chr","", df_mapped$chr_hg19),
    pos_hg19      = df_mapped$pos_hg19,
    rsids         = ifelse(df_mapped$rsids %in% c("", "NA"), NA, df_mapped$rsids)
  ) %>%
    arrange(se) %>% distinct(SNP, .keep_all = TRUE)  # keep best if duplicate SNPs
}

# =========================
# CHUNK 2 — Run liftover + preFormat on all input files
# =========================
in_files <- list.files(IN_DIR, pattern = "_CLEANED\\.txt$", full.names = TRUE)
log_rows <- list()

for (f in in_files) {
  message("Processing: ", basename(f))
  std38 <- standardize_ukb_v38(f)
  lo    <- liftover_to_hg19(std38)
  
  # write hg19 raw with chr/pos
  raw_out <- file.path(OUT_DIR, sub("_CLEANED\\.txt$", "_hg19_raw.tsv", basename(f)))
  fwrite(lo$mapped, raw_out, sep = "\t")
  
  # write preFormat
  pre  <- build_preformat(lo$mapped)
  pre_out <- file.path(OUT_DIR, sub("_CLEANED\\.txt$", "_hg19_preFormat.tsv", basename(f)))
  fwrite(pre, pre_out, sep = "\t")
  
  log_rows[[length(log_rows)+1]] <- data.frame(
    file = basename(f),
    out_raw_hg19 = basename(raw_out),
    out_preFormat_hg19 = basename(pre_out),
    total = lo$summary$total,
    mapped_1to1 = lo$summary$onetoone,
    multi = lo$summary$multi,
    unmapped = lo$summary$unmapped,
    success_pct = round(100*lo$summary$onetoone/lo$summary$total, 2),
    stringsAsFactors = FALSE
  )
}
if (length(log_rows)) {
  liftover_summary <- bind_rows(log_rows)
  fwrite(liftover_summary, file.path(OUT_DIR, "liftover_summary_UKB.csv"))
  print(liftover_summary)
}

# =========================
# CHUNK 3 — cis-window extraction per gene (±500 kb)
# =========================
cis_dir <- file.path(OUT_DIR, "cis500kb_filtered")
dir.create(cis_dir, showWarnings = FALSE)

# hg19 (GRCh37) coordinates — consistent with your Ferkingstad block
genes <- tribble(
  ~gene,   ~chr, ~start,      ~end,
  "IL1RA",    2, 113864791, 113891593,
  "IL1B",     2, 113587328, 113594480,
  "IL6",      7,  22765503,  22771621,
  "IL6R",     1, 154377669, 154441926,
  "TNF",      6,  31543344,  31546113,
  "TNFR1",   12,   6437923,   6451280,
  "TNFR2",    1,  12227060,  12269285
  # "PLAUR",  19,  44150247,  44174699
)

# find a preFormat file for a gene (case-insensitive)
find_preformat_for_gene <- function(g) {
  pat <- paste0("^", g, ".*_hg19_preFormat\\.tsv$")
  hits <- list.files(OUT_DIR, pattern = pat, full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) NA_character_ else hits[1]
}

cis_rows <- pmap_dfr(genes, function(gene, chr, start, end) {
  f <- find_preformat_for_gene(gene)
  if (is.na(f)) return(tibble(gene=gene, n_cis=0, cis_file=NA_character_))
  d <- fread(f)
  d$chr_hg19 <- as.character(d$chr_hg19)
  d$pos_hg19 <- as.integer(d$pos_hg19)
  
  w <- 500*1000L
  s <- max(1L, start - w)
  e <- end + w
  
  cis <- d %>%
    filter(chr_hg19 == as.character(chr),
           !is.na(pos_hg19),
           pos_hg19 >= s, pos_hg19 <= e)
  
  out_path <- file.path(cis_dir, paste0(gene, "_cis500kb_filtered.tsv"))
  fwrite(cis, out_path, sep = "\t")
  tibble(gene=gene, n_cis=nrow(cis), cis_file=out_path)
})

cis_rows
#> per-gene cis tables written to OUT_DIR/cis500kb_filtered

# =========================
# CHUNK 4 — TwoSampleMR exposure format (from cis tables; UNCLUMPED)
# =========================
tsmr_dir <- file.path(OUT_DIR, "TwoSampleMR_exposure")
dir.create(tsmr_dir, showWarnings = FALSE)

exp_rows <- pmap_dfr(cis_rows, function(gene, n_cis, cis_file) {
  if (is.na(cis_file) || !file.exists(cis_file) || n_cis == 0)
    return(tibble(gene=gene, exposure_file=NA_character_, n_kept=0))
  
  d <- fread(cis_file, data.table = FALSE, check.names = FALSE)
  need <- c("SNP","beta","se","effect_allele","other_allele","pval","samplesize")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    warning("Missing columns in ", basename(cis_file), ": ", paste(miss, collapse=", "))
    return(tibble(gene=gene, exposure_file=NA_character_, n_kept=0))
  }
  
  # uppercase alleles
  d$effect_allele <- toupper(d$effect_allele)
  d$other_allele  <- toupper(d$other_allele)
  
  exp <- TwoSampleMR::format_data(
    d %>% transmute(SNP, beta, se,
                    effect_allele, other_allele,
                    eaf = if ("eaf" %in% names(d)) eaf else NA_real_,
                    pval, samplesize),
    type = "exposure",
    snp_col="SNP", beta_col="beta", se_col="se",
    effect_allele_col="effect_allele", other_allele_col="other_allele",
    eaf_col="eaf", pval_col="pval", samplesize_col="samplesize"
  )
  exp$exposure <- gene; exp$id.exposure <- gene
  
  # deduplicate (keep most significant)
  if ("pval.exposure" %in% names(exp)) {
    exp <- exp %>% arrange(pval.exposure) %>% distinct(SNP, .keep_all = TRUE)
  } else {
    exp <- exp %>% distinct(SNP, .keep_all = TRUE)
  }
  
  out <- file.path(tsmr_dir, paste0(gene, "_TwoSampleMR_ready_unclumped.tsv"))
  fwrite(exp, out, sep = "\t")
  tibble(gene=gene, exposure_file=out, n_kept=nrow(exp))
})

exp_rows
#> TwoSampleMR-ready exposure files written to OUT_DIR/TwoSampleMR_exposure

# =========================
# CHUNK 5 — Filter + LD clump cis exposures (p < 1e-5; r2<=0.1; 10,000 kb)
# =========================
clump_dir <- file.path(OUT_DIR, "Exposure_Formatted_Cis")
dir.create(clump_dir, showWarnings = FALSE)

cis_files <- list.files(cis_dir, pattern = "_cis500kb_filtered\\.tsv$", full.names = TRUE)

for (file in cis_files) {
  trait <- sub("_cis500kb_filtered\\.tsv$", "", basename(file))
  d <- fread(file, data.table = FALSE, check.names = FALSE)
  
  # keep biallelic SNPs
  d$effect_allele <- toupper(d$effect_allele)
  d$other_allele  <- toupper(d$other_allele)
  d <- d %>% filter(nchar(effect_allele)==1, nchar(other_allele)==1,
                    effect_allele %in% c("A","C","G","T"),
                    other_allele  %in% c("A","C","G","T"))
  
  d$pval <- to_num(d$pval)
  d <- d %>% filter(!is.na(pval) & pval < 1e-5)
  if (!nrow(d)) { message(" No SNPs passed p < 1e-5 for: ", trait); next }
  
  if (!"rsid" %in% names(d)) d$rsid <- d$SNP
  
  clumped <- tryCatch({
    ieugwasr::ld_clump(
      d = d,
      clump_kb = 10000,
      clump_r2 = 0.1,
      clump_p = 1,
      bfile = REF_BFILE,
      plink_bin = PLINK_BIN
    )
  }, error = function(e) { message(" Clumping failed for ", trait, ": ", e$message); NULL })
  
  if (is.null(clumped) || !nrow(clumped)) { message(" Clumping returned no SNPs for: ", trait); next }
  
  # format for MR (exposure)
  expc <- TwoSampleMR::format_data(
    clumped,
    type = "exposure",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize"
  )
  expc$exposure <- trait; expc$id.exposure <- trait
  
  outc <- file.path(clump_dir, paste0(trait, "_TwoSampleMR_ready_clumped.txt"))
  fwrite(expc, outc, sep = "\t")
  message(" Clumped and saved: ", trait, " | SNPs: ", nrow(expc))
}

# =========================
# CHUNK 6 — COJO formatting (optional)
# =========================
cojo_dir <- file.path(OUT_DIR, "COJO_ready")
dir.create(cojo_dir, showWarnings = FALSE)

clumped_files <- list.files(clump_dir, pattern = "_TwoSampleMR_ready_clumped\\.txt$", full.names = TRUE)

for (file in clumped_files) {
  trait <- sub("_TwoSampleMR_ready_clumped\\.txt$", "", basename(file))
  df <- fread(file, data.table = FALSE)
  
  need <- c("SNP","effect_allele.exposure","other_allele.exposure",
            "eaf.exposure","beta.exposure","se.exposure","pval.exposure","samplesize.exposure")
  if (!all(need %in% names(df))) { warning("Skipping ", trait, " (missing cols)"); next }
  
  cojo <- df %>%
    transmute(
      SNP  = SNP,
      A1   = effect_allele.exposure,
      A2   = other_allele.exposure,
      freq = ifelse(is.na(eaf.exposure), 0.01, eaf.exposure),  # dummy if NA
      b    = beta.exposure,
      se   = se.exposure,
      p    = pval.exposure,
      N    = samplesize.exposure
    )
  
  outfile <- file.path(cojo_dir, paste0(trait, "_COJO_ready.ma"))
  fwrite(cojo, outfile, sep = "\t")
  message(" Saved COJO file: ", basename(outfile))
}

# =========================
# CHUNK 7 — (optional) Merge COJO-selected leads with clumped instruments
# =========================
# If you later produce *_slct.jma.cojo files, you can add the same merge
# routine from Ferkingstad here, pointing to OUT_DIR paths.
