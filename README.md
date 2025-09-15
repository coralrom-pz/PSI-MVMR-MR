PSI-MVMR-MR
Reproducible pipelines for PSI, MR, MVMR, COJO, and COLOC used in our study, “Inflammation mediates the association between social isolation and mental illness: a Mendelian randomization study.”

# 1) Clone
git clone https://github.com/coralrom-pz/PSI-MVMR-MR.git
cd PSI-MVMR-MR

# 2) Recreate the exact R environment
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org"); renv::restore()'

# 3) Repository layout

PSI-MVMR-MR/

├── configs/                       # Repository configuration

├── scripts/                       # R / Rmd / shell pipelines

│   ├── cojo/                      # COJO region builders + post-processing

│   ├── coloc/                     # COLOC per-gene runner + sensitivity + combiner

│   ├── mvmr/                      # MVMR function + batch runner

│   ├── mr/                        # Two-sample MR function + batch runner

│   └── harmonization/             # Exposure/outcome cleaning + formatting

├── bin/                           # Helper CLI tools

│   └── ukbppp_add_rsid.sh         # Adds RSIDs via dbSNP VCF (bcftools)

├── data/                          # Large inputs (NOT versioned)

│   ├── eqtl/

│   │   ├── cis/{Exposure formatted, Outcome formatted, Raw files}

│   │   ├── trans/{Exposure formatted, Outcome formatted, Raw files}

│   │   └── Instrumental Variables_CIS/

│   ├── gwas/

│   │   ├── clumped/{SMI IVs, PSI IVs}

│   │   ├── {Exposure formatted, Outcome formatted}

│   │   └── Raw files/{PSI, SCZ, MDD, etc.}

│   ├── pqtl/

│   │   ├── Ferkingstad/{Exposure formatted, Outcome formatted, Instrumental Variables, Raw files_GRCh38}

│   │   ├── Folkersen/{Exposure formatted, Outcome formatted, Instrumental Variables, RAW files_GWAS catalog}

│   │   └── UKB-PPP/{Exposure formatted, Outcome formatted, Instrumental Variables, Raw Hg19 Files}

│   ├── EU 1000genomes MAF/

│   ├── LD Clumping Data/

│   └── dbSNP/

├── renv.lock                      # Frozen R packages for reproducibility

├── renv/activate.R                # renv bootstrap

└── README.md


# 4) Reproducibility
R packages are locked in renv.lock. Run renv::restore() to recreate the exact environment.

Scripts avoid touching global options; all paths are explicit and commented.

# 5) System dependencies (macOS)
**core tools for dbSNP annotator**
brew install bcftools htslib

**for devtools/gert (optional for development)**
brew install libgit2

**for textshaping/ragg (graphics stack)**
brew install harfbuzz fribidi freetype libpng

# 6) License
GNU General Public License (GPL) v3

# 7) Contact
Questions or issues? Open a GitHub issue or reach out to Coral Romualdo-Perez.
