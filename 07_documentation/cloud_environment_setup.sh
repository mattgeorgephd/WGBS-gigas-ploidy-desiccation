#!/bin/bash
# R toolchain for the WGBS-gigas-ploidy-desiccation pipeline: R 4.3 (Ubuntu 24.04) + Bioconductor 3.18.
# System-wide installs only; nothing here touches a repository checkout. Safe to re-run.
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq || echo "WARN: apt-get update failed" >&2
apt-get install -y --no-install-recommends r-base-core r-base-dev \
  libcurl4-openssl-dev libbz2-dev liblzma-dev zlib1g-dev libglpk40 libxml2-dev libssl-dev \
  r-bioc-genomicranges r-bioc-rtracklayer r-bioc-go.db r-bioc-annotationdbi r-bioc-qvalue \
  r-bioc-rsamtools r-bioc-rhtslib r-bioc-limma r-bioc-complexheatmap r-bioc-gosemsim \
  r-cran-biocmanager r-cran-tidyverse r-cran-data.table r-cran-readxl r-cran-pheatmap \
  r-cran-vegan r-cran-here r-cran-multcomp r-cran-emmeans r-cran-patchwork r-cran-car \
  r-cran-r.utils r-cran-ggrepel r-cran-renv r-cran-knitr r-cran-rmarkdown r-cran-writexl \
  r-cran-openxlsx r-cran-mclust r-cran-emdbook r-cran-gtools r-cran-rcpp r-cran-testthat \
  r-cran-digest r-cran-curl r-cran-optparse r-cran-yaml r-cran-jsonlite r-cran-gplots \
  r-cran-rcolorbrewer r-cran-lme4 r-cran-cowplot r-cran-scales r-cran-stringr \
  || echo "WARN: apt install of R failed" >&2
# Packages Ubuntu does not ship: Bioconductor 3.18 / CRAN (Posit PPM binaries where available).
cat > /tmp/wgbs_r_packages.R <<'EOF'
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"), Ncpus = 4, timeout = 600)
ua <- paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), ua))
lib <- "/usr/local/lib/R/site-library"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
pk <- c("methylKit", "topGO", "simplifyEnrichment", "rrvgo", "agricolae", "Biostrings", "ggbeeswarm")
need <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) BiocManager::install(need, lib = lib, update = FALSE, ask = FALSE)
EOF
Rscript /tmp/wgbs_r_packages.R || echo "WARN: R package install failed" >&2
