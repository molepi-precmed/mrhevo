## compute_alpha_ukb.R
## Compute instrument-exposure coefficient alpha_hat for each locus
## from UK Biobank Olink proteomics GWAS summary statistics for PDCD1.
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
##
## Bim used for: (chr, pos hg38) -> rsid lookup; allele alignment done vs 1000G LD A1.

library(data.table)
library(jsonlite)
source("ld_functions.R")   # for define_loci_expanded only

## ---- Parameters ----
tar_file       <- "pdcd1_stats/PDCD1_Q15116_OID21396_v1_Oncology.tar"
bim_dir        <- "refpop/bim_by_chr"
## ---- LD reference panel paths (on genoscores) ----
zarr_dir_ukb <- "/opt/datastore/genome/LD_Eur_18mvariants/int8"
zarr_dir_1kg <- "/opt/datastore/genome/1000G/1kg_ld_eur/ld_1kg_eur"
zarr_dir     <- zarr_dir_ukb   # UKBB EUR zarr (362k samples)
genoscores_host <- "pmckeigue@genoscores.cphs.mvm.ed.ac.uk"
gap_mb            <- 1.0
window_mb         <- 1.0
min_eig_frac      <- 0.01
p_hit             <- 1e-6
p_cand            <- 1e-5
exclude_hla       <- TRUE
info_threshold    <- 0.3
min_maf_threshold <- 0.02
out_file       <- "alpha_pdcd1_UKB.rds"

## ---- Profiling helpers ----
.t_script <- proc.time()["elapsed"]
.checkpoints <- list()
checkpoint <- function(label) {
    elapsed <- proc.time()["elapsed"] - .t_script
    mem_mb  <- sum(gc(verbose = FALSE)[, 2L])
    .checkpoints[[label]] <<- list(elapsed_s = elapsed, mem_mb = mem_mb)
}
checkpoint("start")

## ---- Step 1: Load summary statistics ----
message("Loading summary statistics from tar archive...")
prefilter_log10p <- 5   # p < 1e-5 prefilter during tar read

members <- system2("tar", c("-tf", tar_file), stdout = TRUE)
members <- members[grepl("\\.gz$", members)]

read_one <- function(m) {
    cmd <- sprintf("tar -xOf '%s' '%s' | zcat", tar_file, m)
    dt  <- fread(cmd = cmd,
                 select = c("CHROM", "GENPOS", "ALLELE0", "ALLELE1",
                            "A1FREQ", "INFO", "N", "BETA", "SE", "LOG10P"),
                 colClasses = list(character = c("ALLELE0", "ALLELE1")))
    dt[LOG10P > prefilter_log10p]
}
stats <- rbindlist(lapply(members, read_one))

setnames(stats,
         c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ",
           "N",               "BETA",    "SE",      "LOG10P"),
         c("chr",   "pos",    "Allele2", "Allele1", "Freq1",
           "TotalSampleSize", "Effect",  "StdErr",  "log10p"))
stats[, Allele1 := toupper(Allele1)]
stats[, Allele2 := toupper(Allele2)]
stats[, chr := as.integer(chr)]
setorder(stats, chr, pos)
message(sprintf("  %d SNPs pass pre-filter p < 1e-5 (LOG10P > 5)", nrow(stats)))
checkpoint("after_load_stats")

## ---- Step 1b: Filter SNPs ----
message("Filtering SNPs...")
stats <- stats[INFO >= info_threshold]
stats[, INFO := NULL]
message(sprintf("  %d SNPs pass INFO >= %.1f", nrow(stats), info_threshold))

## Convert log10p to Pvalue for define_loci_expanded
stats[, Pvalue := 10^(-log10p)]
stats[, log10p := NULL]

is_ambiguous <- function(x, y)
    (x == "A" & y == "T") | (x == "T" & y == "A") |
    (x == "C" & y == "G") | (x == "G" & y == "C")
n_before <- nrow(stats)
stats <- stats[!is_ambiguous(Allele1, Allele2)]
message(sprintf("  %d SNPs after dropping strand-ambiguous (%d dropped)",
                nrow(stats), n_before - nrow(stats)))
stats <- stats[pmin(Freq1, 1 - Freq1) >= min_maf_threshold]
message(sprintf("  %d SNPs after MAF >= %.2f filter", nrow(stats), min_maf_threshold))
checkpoint("after_filter")

## ---- Step 2: Define loci (Zhou et al. locus expansion) ----
message("Defining loci (hit-anchored expansion)...")
stats <- define_loci_expanded(stats, p_hit = p_hit, p_cand = p_cand,
                               window_mb = window_mb, gap_mb = gap_mb,
                               exclude_hla = exclude_hla)
message(sprintf("  %d SNPs in %d loci", nrow(stats), uniqueN(stats$locus_id)))

## Save top SNP per locus for gene lookup
top_snps_ukb <- stats[stats[, .I[which.min(Pvalue)], by = locus_id]$V1,
                      .(locus_id, chr, top_pos = pos)]
stats[, Pvalue := NULL]
checkpoint("after_define_loci")

## ---- Step 3: (chr, pos hg38) -> rsid via per-chromosome bim RDS files ----
## Bim used solely as coordinate lookup; allele alignment done vs 1000G LD below.
message("Loading per-chr bim files for (chr, pos) -> rsid lookup...")
bim_chrs <- unique(stats$chr)
stats <- rbindlist(lapply(bim_chrs, function(ch) {
    f <- file.path(bim_dir, sprintf("chr_%s.rds", ch))
    if (!file.exists(f)) return(NULL)
    bim_ch <- readRDS(f)[, .(chr, rsid, pos)]   # one chr at a time; freed after join
    bim_ch[stats[chr == ch], on = c("chr", "pos")]
}))
checkpoint("after_bim_load")

n_before <- nrow(stats)
stats <- stats[!is.na(rsid)]
stats <- unique(stats, by = c("chr", "pos"))
message(sprintf("  %d SNPs with rsid from bim (%d dropped: no bim match or multi-allelic)",
                nrow(stats), n_before - nrow(stats)))

## Recompute loci on matched set
stats[, prev_pos  := shift(pos, 1L, type = "lag"), by = chr]
stats[, new_locus := is.na(prev_pos) | (pos - prev_pos) > gap_mb * 1e6]
stats[, locus_n   := cumsum(new_locus), by = chr]
stats[, locus_id  := paste0("chr", chr, "_locus", locus_n)]
stats[, c("prev_pos", "new_locus", "locus_n") := NULL]
message(sprintf("  %d loci after bim matching", uniqueN(stats$locus_id)))

## Minimum minor allele frequency per locus (from UKBB proteomics summary stats)
locus_min_maf <- stats[, .(min_maf = min(pmin(Freq1, 1 - Freq1))), by = locus_id]

## Remap top_snps_ukb locus_id to post-match locus_id
top_snps_ukb <- merge(top_snps_ukb,
                      stats[, .(chr, pos, locus_id_matched = locus_id)],
                      by.x = c("chr", "top_pos"), by.y = c("chr", "pos"),
                      all.x = FALSE)
top_snps_ukb[, locus_id := locus_id_matched][, locus_id_matched := NULL]

checkpoint("after_bim_match")

## ---- Step 4: Serialize per-locus SNP data and run zarr LD on genoscores ----
message("Building per-locus JSON for zarr LD computation on genoscores...")
locus_ids <- unique(stats$locus_id)
single_in_stats <- stats[, .N, by = locus_id][N == 1L, locus_id]

job <- list(
    zarr_dir     = zarr_dir,
    min_eig_frac = min_eig_frac,
    loci = lapply(locus_ids, function(lid) {
        snps <- stats[locus_id == lid]
        list(
            locus_id = lid,
            chr      = snps$chr[1L],
            snps     = lapply(seq_len(nrow(snps)), function(i) list(
                rsid            = snps$rsid[i],
                Allele1         = snps$Allele1[i],
                Allele2         = snps$Allele2[i],
                Effect          = snps$Effect[i],
                Freq1           = snps$Freq1[i],
                StdErr          = snps$StdErr[i],
                TotalSampleSize = snps$TotalSampleSize[i],
                pos             = snps$pos[i]
            ))
        )
    })
)

json_in  <- tempfile(fileext = ".json")
json_err <- tempfile(fileext = ".err")
json_out <- tempfile(fileext = ".json")
write_json(job, json_in, auto_unbox = TRUE, digits = 8)
on.exit({ unlink(c(json_in, json_err, json_out)) }, add = TRUE)

message(sprintf("Sending %d loci to genoscores via SSH...", length(locus_ids)))
ret <- system(sprintf(
    "scp -q genoscores_zarr_loci.py '%s':~/ && ssh '%s' 'python3 ~/genoscores_zarr_loci.py' < '%s' > '%s' 2>'%s'",
    genoscores_host, genoscores_host, json_in, json_out, json_err
))
if (file.exists(json_err) && file.size(json_err) > 0L)
    message(paste(readLines(json_err, warn = FALSE), collapse = "\n"))
if (ret != 0L) stop("genoscores_zarr_loci.py failed on genoscores (exit ", ret, ")")
checkpoint("after_loci")

response <- read_json(json_out, simplifyVector = FALSE)$loci
alpha_dt <- rbindlist(lapply(response, function(r) data.table(
    qtlname      = r$locus_id,
    chr          = as.character(r$chr),
    locus_start  = as.integer(r$locus_start),
    locus_end    = as.integer(r$locus_end),
    n_snps       = as.integer(r$n_snps),
    eff_rank     = as.integer(r$eff_rank),
    n_ld_blocks  = as.integer(r$n_ld_blocks),
    alpha_hat    = as.numeric(r$alpha_hat),
    se.alpha_hat = as.numeric(r$se_alpha_hat),
    sd.Z         = as.numeric(r$sd_Z)
)))
setorder(alpha_dt, chr, locus_start)
alpha_dt <- merge(alpha_dt, locus_min_maf, by.x = "qtlname", by.y = "locus_id",
                  all.x = TRUE)

## ---- Step 5: Report and save ----
message(sprintf("\nCompleted: %d loci", nrow(alpha_dt)))
print(alpha_dt[, .(qtlname, chr, n_snps, eff_rank, alpha_hat, se.alpha_hat)])

single_snp <- alpha_dt[n_snps == 1L & qtlname %in% single_in_stats]
if (nrow(single_snp) > 0L) {
    single_stats <- stats[locus_id %in% single_snp$qtlname,
                          .(qtlname = locus_id, expected = abs(Effect))]
    check <- merge(single_snp[, .(qtlname, alpha_hat)], single_stats, by = "qtlname")
    max_rel <- check[, max(abs(alpha_hat - expected) / expected)]
    message(sprintf("Single-SNP check: max |alpha_hat - |Effect|| / |Effect| = %.2e", max_rel))
} else {
    message("Single-SNP check: no single-SNP loci in both stats and LD matrix")
}

saveRDS(alpha_dt, out_file)
message(sprintf("Saved: %s", out_file))

## ---- Step 6: Nearest protein-coding gene to top SNP per locus ----
message("Looking up nearest protein-coding gene per locus...")
genes <- fread("refGene.txt.gz", header = FALSE,
               select = c(2L, 3L, 5L, 6L, 13L),
               col.names = c("accession", "chrom", "txStart", "txEnd", "symbol"))
genes[, chr_g := suppressWarnings(as.integer(sub("^chr", "", chrom)))]
genes <- genes[!is.na(chr_g) & startsWith(accession, "NM_")]
genes <- genes[, .(txStart = min(txStart), txEnd = max(txEnd)),
               by = .(chr = chr_g, symbol)]

top_snps_ukb <- top_snps_ukb[locus_id %in% alpha_dt$qtlname]

find_nearest <- function(chr_val, pos_val) {
    g <- genes[chr == chr_val]
    if (nrow(g) == 0L) return(NA_character_)
    dist <- pmax(0L, pmax(g$txStart - pos_val, pos_val - g$txEnd))
    g$symbol[which.min(dist)]
}
top_snps_ukb[, gene := mapply(find_nearest, chr, top_pos)]

locus_genes_ukb <- top_snps_ukb[, .(qtlname = locus_id, gene)]
saveRDS(locus_genes_ukb, "locus_genes_ukb.rds")
message("Saved: locus_genes_ukb.rds")

## ---- Profiling summary ----
cp_names <- names(.checkpoints)
cp_dt <- data.table(
    step      = cp_names,
    elapsed_s = sapply(.checkpoints, `[[`, "elapsed_s"),
    mem_mb    = sapply(.checkpoints, `[[`, "mem_mb")
)
cp_dt[, step_s := c(elapsed_s[1L], diff(elapsed_s))]

message("\n=== Step timing and memory ===")
print(cp_dt[, .(step, step_s = round(step_s, 1), cumul_s = round(elapsed_s, 1),
                mem_mb = round(mem_mb, 0))])

message(sprintf("  SSH + zarr LD total: %.1fs for %d loci",
                .checkpoints[["after_loci"]]$elapsed_s - .checkpoints[["after_bim_match"]]$elapsed_s,
                length(locus_ids)))
