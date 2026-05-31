## compute_gamma.R
## Compute instrument-outcome coefficient gamma_hat for each locus
## from GWAS summary statistics (z-scores) for SNP-T1D associations.
##
## LD source: 1000G EUR panel (refpop/kg.2020.hg38.eur, 503 samples, hg38).
## Locus definition: Zhou et al. procedure — hits at p < 1e-6, candidates
##   at p < 1e-5 within 1 Mb of any hit (see define_loci_expanded in ld_functions.R).
##
## Method: theorymethods.Rmd
##   gamma_hat = alpha_hat * (gamma_u*^T alpha_m) / Var(S)   [eq. 16]
##   SE = 1 / sqrt(N * Var(Z) * p*(1-p))                     [eq. 17]
##
## Source of T1D GWAS: Iakovliev et al. 2023 AJHG, PMID 37164005
##
## Bim used for allele alignment of PDCD1 and T1D z-scores to 1000G a1 convention.
## Both allele alignment and LD A1 are from the same 1000G bim, so the second
## alignment step (bim_a1 -> LD_A1) is always a no-op, but handled correctly.

library(data.table)
library(BEDMatrix)
source("ld_functions.R")

## ---- Parameters ----
stats_file   <- "pdcd1_stats/OID00791_OID21396_SCALLOP_UKBB_MA_rsidannotated_filtered.tsv"
t1d_file     <- "t1d_stats/data/gwasresults.RData.gz"
bim_dir        <- "refpop/bim_by_chr"
bed_file       <- "refpop/kg.2020.hg38.eur"
gap_mb         <- 1.0
window_mb      <- 1.0
min_eig_frac   <- 0.01
p_hit          <- 1e-6
p_cand         <- 1e-5
exclude_hla    <- TRUE
out_file       <- "gamma_t1d.rds"

## T1D GWAS sample size (Iakovliev et al. 2023)
N_cases  <- 4964L
N_ctrls  <- 7497L
N_gamma  <- N_cases + N_ctrls
p_case   <- N_cases / N_gamma

## ---- Step 1: Load T1D GWAS summary statistics ----
## z-scores assumed coded for bim a1 allele (plink convention).
message("Loading T1D GWAS summary statistics...")
tmp_rdata <- tempfile(fileext = ".RData")
system2("gunzip", args = c("-c", t1d_file), stdout = tmp_rdata)
load(tmp_rdata)
unlink(tmp_rdata)
t1d <- as.data.table(results); rm(results)
t1d[, CHR := as.integer(as.character(CHR))]
setnames(t1d, c("SNP", "CHR", "z"), c("rsid", "chr", "z_t1d"))
t1d[, c("BP", "minuslog10p") := NULL]

if (anyDuplicated(t1d$rsid))
    t1d <- t1d[t1d[, .I[which.max(abs(z_t1d))], by = rsid]$V1]
setkey(t1d, rsid)
message(sprintf("  %d T1D SNPs loaded", nrow(t1d)))

## ---- Step 2: Load PDCD1 summary statistics ----
message("Loading PDCD1 summary statistics...")
stats <- fread(stats_file,
               select = c("chr", "pos", "rsid", "Allele1", "Allele2",
                          "Freq1", "Effect", "StdErr", "TotalSampleSize",
                          "Pvalue"))
stats[, Allele1 := toupper(Allele1)]
stats[, Allele2 := toupper(Allele2)]
setorder(stats, chr, pos)
message(sprintf("  %d PDCD1 SNPs loaded", nrow(stats)))

## Drop strand-ambiguous SNPs
is_ambiguous <- function(x, y)
    (x == "A" & y == "T") | (x == "T" & y == "A") |
    (x == "C" & y == "G") | (x == "G" & y == "C")
stats <- stats[!is_ambiguous(Allele1, Allele2)]

## ---- Step 3: Locus expansion (Zhou et al.) ----
message("Defining loci (hit-anchored expansion)...")
stats <- define_loci_expanded(stats, p_hit = p_hit, p_cand = p_cand,
                               window_mb = window_mb, gap_mb = gap_mb,
                               exclude_hla = exclude_hla)
message(sprintf("  %d PDCD1 SNPs in %d loci", nrow(stats), uniqueN(stats$locus_id)))

## ---- Step 4: Align PDCD1 and T1D to 1000G bim a1 convention ----
## Bim a1 convention is the same as 1000G LD A1, so this single alignment step
## ensures both PDCD1 Effect and T1D z are coded for the same allele before LD lookup.
message("Loading per-chr bim files for allele alignment...")
bim_chrs <- unique(stats$chr)
bim <- rbindlist(lapply(bim_chrs, function(ch) {
    f <- file.path(bim_dir, sprintf("chr_%s.rds", ch))
    if (!file.exists(f)) return(NULL)
    readRDS(f)
}))
setkey(bim, rsid)

## Align PDCD1 to bim a1
stats <- bim[, .(rsid, bim_a1 = a1, bim_a2 = a2)][stats, on = "rsid"]
stats[, allele_match := fcase(
    Allele1 == bim_a1 & Allele2 == bim_a2, "ok",
    Allele1 == bim_a2 & Allele2 == bim_a1, "flip",
    default = "drop"
)]
stats[allele_match == "flip", `:=`(Effect = -Effect, Freq1 = 1 - Freq1,
                                    Allele1 = bim_a1, Allele2 = bim_a2)]
n_before <- nrow(stats)
stats <- stats[!is.na(bim_a1) & allele_match != "drop"]
message(sprintf("  %d PDCD1 SNPs aligned to bim (%d dropped)", nrow(stats), n_before - nrow(stats)))
## After this: Allele1 == bim_a1 == LD A1 for all remaining stats rows.
## T1D z-scores are coded for bim_a1 — consistent with stats$Allele1.

locus_n_alpha <- stats[, .(n_snps_alpha = .N), by = locus_id]

## ---- Step 5: Open 1000G BED and process each locus ----
message("Opening 1000G BED file (memory-mapped)...")
bed    <- BEDMatrix(bed_file)
get_ld <- function(chrom, rsids) get_locus_ld_1kg(bed, bim_dir, chrom, rsids)

message("Processing loci...")
locus_ids <- unique(stats$locus_id)

process_locus <- function(lid) {
    alpha_snps <- stats[locus_id == lid]

    ## Restrict to rsids present in T1D GWAS
    candidate_rsids <- intersect(alpha_snps$rsid, t1d$rsid)
    if (length(candidate_rsids) == 0L) return(NULL)
    alpha_snps <- alpha_snps[rsid %in% candidate_rsids]

    ## Build named z-score vector (names = rsid) for locus_compute_gamma
    z_gamma <- setNames(t1d[.(candidate_rsids)]$z_t1d, candidate_rsids)

    res <- locus_compute_gamma(
        chrom        = alpha_snps$chr[1L],
        snps         = alpha_snps,
        z_gamma      = z_gamma,
        N_gamma      = N_gamma,
        p_case       = p_case,
        get_ld_fn    = get_ld,
        min_eig_frac = min_eig_frac
    )
    if (is.null(res)) return(NULL)
    if (res$eff_rank < res$n_snps / 10)
        message(sprintf("  %s: n_snps=%d, eff_rank=%d (high LD, truncated)",
                        lid, res$n_snps, res$eff_rank))
    dt <- as.data.table(c(list(qtlname = lid), res))
    setnames(dt, "n_snps", "n_snps_gamma")
    dt
}

results_list <- lapply(locus_ids, process_locus)
gamma_dt     <- rbindlist(Filter(Negate(is.null), results_list))
setorder(gamma_dt, chr, locus_start)

gamma_dt <- locus_n_alpha[gamma_dt, on = .(locus_id = qtlname)]
setnames(gamma_dt, "locus_id", "qtlname")
setcolorder(gamma_dt, c("qtlname", "chr", "locus_start", "locus_end",
                         "n_snps_alpha", "n_snps_gamma", "eff_rank",
                         "gamma_hat", "se.gamma_hat"))

## ---- Step 6: Report and save ----
message(sprintf("\nCompleted: %d loci with gamma estimates (of %d PDCD1 loci)",
                nrow(gamma_dt), length(locus_ids)))
n_partial <- gamma_dt[n_snps_gamma < n_snps_alpha, .N]
if (n_partial > 0L)
    message(sprintf("  %d loci have fewer T1D-matched SNPs than PDCD1 SNPs", n_partial))
print(gamma_dt[, .(qtlname, chr, n_snps_alpha, n_snps_gamma, eff_rank,
                    gamma_hat, se.gamma_hat)])

saveRDS(gamma_dt, out_file)
message(sprintf("Saved: %s", out_file))

## ---- Step 7: Nearest protein-coding gene to top SNP per locus ----
message("Looking up nearest protein-coding gene per locus...")
genes <- fread("refGene.txt.gz", header = FALSE,
               select = c(2L, 3L, 5L, 6L, 13L),
               col.names = c("accession", "chrom", "txStart", "txEnd", "symbol"))
genes[, chr_g := suppressWarnings(as.integer(sub("^chr", "", chrom)))]
genes <- genes[!is.na(chr_g) & startsWith(accession, "NM_")]
genes <- genes[, .(txStart = min(txStart), txEnd = max(txEnd)),
               by = .(chr = chr_g, symbol)]

top_snps <- stats[stats[, .I[which.min(Pvalue)], by = locus_id]$V1,
                  .(locus_id, chr, top_pos = pos)]
top_snps <- top_snps[locus_id %in% gamma_dt$qtlname]

find_nearest <- function(chr_val, pos_val) {
    g <- genes[chr == chr_val]
    if (nrow(g) == 0L) return(NA_character_)
    dist <- pmax(0L, pmax(g$txStart - pos_val, pos_val - g$txEnd))
    g$symbol[which.min(dist)]
}
top_snps[, gene := mapply(find_nearest, chr, top_pos)]

locus_genes <- top_snps[, .(qtlname = locus_id, gene)]
saveRDS(locus_genes, "locus_genes.rds")
message("Saved: locus_genes.rds")
