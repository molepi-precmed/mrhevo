## compute_alpha.R
## Compute instrument-exposure coefficient alpha_hat for each locus
## from SCALLOP/UKBB meta-analysis summary statistics for PDCD1.
##
## LD source: 1000G EUR panel (refpop/kg.2020.hg38.eur, 503 samples, hg38).
## Locus definition: Zhou et al. procedure — hits at p < 1e-6, candidates
##   at p < 1e-5 within 1 Mb of any hit (see define_loci_expanded in ld_functions.R).
##
## Method: theorymethods.Rmd section "Constructing scalar instruments from multiple SNPs"
##   alpha_m = C_g^{-1} alpha_u*   (multivariable per-SD coefficients)
##   alpha_hat = ||beta_m||          (per-allele norm, eq. 7)
##   Var(S) = alpha_m^T C_g alpha_m  (analytical)
##   SE via Fisher information (eq. 10)

library(data.table)
library(BEDMatrix)
source("ld_functions.R")

## ---- Parameters ----
stats_file   <- "pdcd1_stats/OID00791_OID21396_SCALLOP_UKBB_MA_rsidannotated_filtered.tsv"
bim_dir        <- "refpop/bim_by_chr"
bed_file       <- "refpop/kg.2020.hg38.eur"
gap_mb         <- 1.0
window_mb      <- 1.0
min_eig_frac   <- 0.01
p_hit          <- 1e-6
p_cand         <- 1e-5
exclude_hla    <- TRUE
out_file       <- "alpha_pdcd1_OID00791.rds"

## ---- Step 1: Load summary statistics ----
message("Loading summary statistics...")
stats <- fread(stats_file,
               select = c("chr", "pos", "rsid", "Allele1", "Allele2",
                          "Freq1", "Effect", "StdErr", "TotalSampleSize",
                          "Pvalue"))
stats[, Allele1 := toupper(Allele1)]
stats[, Allele2 := toupper(Allele2)]
setorder(stats, chr, pos)
message(sprintf("  %d SNPs loaded", nrow(stats)))

## ---- Step 1b: Drop strand-ambiguous SNPs ----
is_ambiguous <- function(x, y)
    (x == "A" & y == "T") | (x == "T" & y == "A") |
    (x == "C" & y == "G") | (x == "G" & y == "C")
n_before <- nrow(stats)
stats <- stats[!is_ambiguous(Allele1, Allele2)]
message(sprintf("  %d SNPs after dropping strand-ambiguous (%d dropped)",
                nrow(stats), n_before - nrow(stats)))

## ---- Step 2: Define loci (Zhou et al. locus expansion) ----
message("Defining loci (hit-anchored expansion)...")
stats <- define_loci_expanded(stats, p_hit = p_hit, p_cand = p_cand,
                               window_mb = window_mb, gap_mb = gap_mb,
                               exclude_hla = exclude_hla)
message(sprintf("  %d SNPs in %d loci (p_hit < %.0e, p_cand < %.0e, window %.1f Mb)",
                nrow(stats), uniqueN(stats$locus_id), p_hit, p_cand, window_mb))

## ---- Step 3: Open 1000G BED (memory-mapped) ----
message("Opening 1000G BED file (memory-mapped)...")
bed    <- BEDMatrix(bed_file)
get_ld <- function(chrom, rsids) get_locus_ld_1kg(bed, bim_dir, chrom, rsids)

## ---- Step 4: Process each locus ----
message("Processing loci...")
locus_ids <- unique(stats$locus_id)

process_locus <- function(lid) {
    snps <- stats[locus_id == lid]
    res  <- locus_compute_alpha(snps$chr[1L], snps, get_ld, min_eig_frac)
    if (is.null(res)) return(NULL)
    if (res$eff_rank < res$n_snps / 10)
        message(sprintf("  %s: n_snps=%d, eff_rank=%d (high LD, truncated)",
                        lid, res$n_snps, res$eff_rank))
    as.data.table(c(list(qtlname = lid), res))
}

results  <- lapply(locus_ids, process_locus)
alpha_dt <- rbindlist(Filter(Negate(is.null), results))
setorder(alpha_dt, chr, locus_start)

## ---- Step 5: Report and save ----
message(sprintf("\nCompleted: %d loci", nrow(alpha_dt)))
print(alpha_dt[, .(qtlname, chr, n_snps, eff_rank, alpha_hat, se.alpha_hat)])

## Sanity check: loci with exactly 1 SNP in both stats and LD: alpha_hat == |Effect|
single_in_stats <- stats[, .N, by = locus_id][N == 1L, locus_id]
single_snp <- alpha_dt[n_snps == 1L & qtlname %in% single_in_stats]
if (nrow(single_snp) > 0L) {
    single_stats <- stats[locus_id %in% single_snp$qtlname,
                          .(qtlname = locus_id, expected = abs(Effect))]
    check <- merge(single_snp[, .(qtlname, alpha_hat)], single_stats, by = "qtlname")
    max_rel <- check[, max(abs(alpha_hat - expected) / expected)]
    message(sprintf("Single-SNP check: max |alpha_hat - |Effect|| / |Effect| = %.2e", max_rel))
} else {
    message("Single-SNP check: no loci with exactly 1 SNP in both stats and LD matrix")
}

saveRDS(alpha_dt, out_file)
message(sprintf("Saved: %s", out_file))
