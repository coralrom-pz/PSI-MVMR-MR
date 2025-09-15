# config/paths.R
#Reproduce environment in terminal with Rscript -e 'renv::restore(prompt = FALSE)'
cfg <- list(
  data = list(
    gwas = list(
      raw      = here("data/gwas/Raw files"),
      clumped  = here("data/gwas/clumped"),
      exp_fmt  = here("data/gwas/Exposure formatted"),
      out_fmt  = here("data/gwas/Outcome formatted")
    ),
    eqtl = list(
      cis = list(
        raw      = here("data/eqtl/cis/Raw files"),
        exp_fmt  = here("data/eqtl/cis/Exposure formatted"),
        out_fmt  = here("data/eqtl/cis/Outcome formatted")
      ),
      trans = list(
        raw      = here("data/eqtl/trans/Raw files"),
        exp_fmt  = here("data/eqtl/trans/Exposure formatted"),
        out_fmt  = here("data/eqtl/trans/Outcome formatted")
      ),
      ivs_cis = here("data/eqtl/Instrumental Variables_CIS")
    ),
    pqtl = list(
      folkersen = list(
        raw      = here("data/pqtl/Folkersen/RAW files_GWAS catalog"),
        exp_fmt  = here("data/pqtl/Folkersen/Exposure formatted"),
        out_fmt  = here("data/pqtl/Folkersen/Outcome formatted"),
        ivs      = here("data/pqtl/Folkersen/Instrumental Variables")
      ),
      ferkingstad = list(
        raw      = here("data/pqtl/Ferkingstad/Raw files_GRCh38"),
        exp_fmt  = here("data/pqtl/Ferkingstad/Exposure formatted"),
        out_fmt  = here("data/pqtl/Ferkingstad/Outcome formatted"),
        ivs      = here("data/pqtl/Ferkingstad/Instrumental Variables")
      ),
      ukbppp = list(
        raw_hg19 = here("data/pqtl/UKB-PPP/Raw Hg19 Files"),
        exp_fmt  = here("data/pqtl/UKB-PPP/Exposure formatted"),
        out_fmt  = here("data/pqtl/UKB-PPP/Outcome formatted"),
        ivs      = here("data/pqtl/UKB-PPP/Instrumental Variables")
      )
    ),
    maf_1kg = here("data/EU 1000genomes MAF"),
    ld_ref  = here("data/LD Clumping Data")
  ),
  outputs = list(
    cojo = list(
      regions_500kb   = here("outputs/cojo/cojo_regions/500kb"),
      regions_500kb_p = here("outputs/cojo/cojo_regions/500kb/pQTLs"),
      results_500kb   = here("outputs/cojo/cojo_results/500kb"),
      cond_results    = here("outputs/cojo/cojo_cond_results"),
      slct_results    = here("outputs/cojo/cojo_slct_results"),
      final_iv        = here("outputs/cojo/final_instruments"),
      final_iv_clean  = here("outputs/cojo/final_instruments_cleaned"),
      combined        = here("outputs/cojo/combined")
    ),
    mvmr  = here("outputs/mvmr"),
    coloc = here("outputs/coloc"),
    figs  = here("outputs/figures")
  ),
  resources = list(
    dbsnp_vcf = here("data/dbSNP/dbsnp.vcf.gz") # keep .tbi next to it
  ),
  constants = list(
    n_gwas_isolation = 452302,
    n_eqtl_default   = 31684,
    n_pqtl_default   = 21758
  )
)
