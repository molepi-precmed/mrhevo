## compute_gamma_ukb.R
## Compute instrument-outcome coefficient gamma_hat for each locus
## using UK Biobank Olink PDCD1 summary stats to define loci,
## and T1D GWAS z-scores (Iakovliev et al. 2023) as the outcome.
##
## LD source: UKBB EUR zarr store on genoscores (362k samples, GRCh37, 18M variants).
##   SNP matching by rsid (build-agnostic); zarr LD accessed via SSH.
## Locus definition: Zhou et al. procedure — hits at p < 1e-6, candidates
##   at p < 1e-5 within 1 Mb of any hit (see define_loci_expanded in ld_functions.R).
##
## Method: theorymethods.Rmd eqs. 16, 17 (analytical gamma_hat and SE).
##
## Bim used for: (1) (chr, pos hg38) -> rsid for PDCD1 SNPs;
##               (2) a1/a2 for allele alignment of PDCD1 and T1D to bim_a1.
##               Bim_a1 re-alignment to zarr LD A1 is handled in run_zarr_loci.py.

library(data.table)
library(jsonlite)
source("ld_functions.R")   # for define_loci_expanded only

## ---- Parameters ----
tar_file       <- "pdcd1_stats/PDCD1_Q15116_OID21396_v1_Oncology.tar"
t1d_file       <- "t1d_stats/data/gwasresults.RData.gz"
bim_dir        <- "refpop/bim_by_chr"
zarr_dir       <- "/opt/datastore/genome/LD_Eur_18mvariants/int8"
genoscores_host <- "pmckeigue@genoscores.cphs.mvm.ed.ac.uk"
gap_mb         <- 1.0
window_mb      <- 1.0
min_eig_frac   <- 0.01
p_hit          <- 1e-6
p_cand         <- 1e-5
exclude_hla    <- TRUE
info_threshold <- 0.3

## T1D GWAS sample size (Iakovliev et al. 2023)
N_cases  <- 4964L
N_ctrls  <- 7497L
N_gamma  <- N_cases + N_ctrls
p_case   <- N_cases / N_gamma

## ---- Profiling helpers ----
.t_script <- proc.time()["elapsed"]
.checkpoints <- list()
checkpoint <- function(label) {
    elapsed <- proc.time()["elapsed"] - .t_script
    mem_mb  <- sum(gc(verbose = FALSE)[, 2L])
    .checkpoints[[label]] <<- list(elapsed_s = elapsed, mem_mb = mem_mb)
}
checkpoint("start")

## ---- Step 1: Load T1D GWAS summary statistics ----
message("Loading T1D GWAS summary statistics...")
tmp_rdata <- tempfile(fileext = ".RData")
system2("gunzip", args = c("-c", t1d_file), stdout = tmp_rdata)
load(tmp_rdata); unlink(tmp_rdata)
t1d <- as.data.table(results); rm(results)
t1d[, CHR := as.integer(as.character(CHR))]
setnames(t1d, c("SNP", "CHR", "z"), c("rsid", "chr", "z_t1d"))
t1d[, c("BP", "minuslog10p") := NULL]
if (anyDuplicated(t1d$rsid))
    t1d <- t1d[t1d[, .I[which.max(abs(z_t1d))], by = rsid]$V1]
setkey(t1d, rsid)
message(sprintf("  %d T1D SNPs loaded", nrow(t1d)))
checkpoint("after_load_t1d")

## ---- Step 2: Load UKB PDCD1 summary statistics ----
message("Loading UKB PDCD1 summary statistics from tar archive...")
prefilter_log10p <- 5   # p < 1e-5 prefilter during tar read
members <- system2("tar", c("-tf", tar_file), stdout = TRUE)
members <- members[grepl("\\.gz$", members)]

read_one <- function(m) {
    dt <- fread(cmd = sprintf("tar -xOf '%s' '%s' | zcat", tar_file, m),
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

stats <- stats[INFO >= info_threshold]; stats[, INFO := NULL]
message(sprintf("  %d UKB PDCD1 SNPs pass INFO >= %.1f and p < 1e-5", nrow(stats), info_threshold))

## Convert log10p to Pvalue for define_loci_expanded
stats[, Pvalue := 10^(-log10p)]
stats[, log10p := NULL]

## Drop strand-ambiguous SNPs
is_ambiguous <- function(x, y)
    (x=="A"&y=="T")|(x=="T"&y=="A")|(x=="C"&y=="G")|(x=="G"&y=="C")
stats <- stats[!is_ambiguous(Allele1, Allele2)]
checkpoint("after_load_pdcd1")

## ---- Step 3: Define loci (Zhou et al. locus expansion) ----
message("Defining loci (hit-anchored expansion)...")
stats <- define_loci_expanded(stats, p_hit = p_hit, p_cand = p_cand,
                               window_mb = window_mb, gap_mb = gap_mb,
                               exclude_hla = exclude_hla)
message(sprintf("  %d UKB PDCD1 SNPs in %d loci", nrow(stats), uniqueN(stats$locus_id)))

## Save top SNP per locus for gene lookup
top_snps_ukb <- stats[stats[, .I[which.min(Pvalue)], by = locus_id]$V1,
                      .(locus_id, chr, top_pos = pos)]
stats[, Pvalue := NULL]

## ---- Step 4: Bim match: (chr, pos) -> rsid and allele alignment ----
message("Loading per-chr bim files for (chr, pos) -> rsid and allele alignment...")
bim_chrs <- unique(stats$chr)
stats <- rbindlist(lapply(bim_chrs, function(ch) {
    f <- file.path(bim_dir, sprintf("chr_%s.rds", ch))
    if (!file.exists(f)) return(NULL)
    bim_ch <- readRDS(f)[, .(chr, pos, rsid, bim_a1 = a1, bim_a2 = a2)]
    bim_ch[stats[chr == ch], on = c("chr", "pos")]
}))
checkpoint("after_bim_load")
stats[, allele_match := fcase(
    Allele1 == bim_a1 & Allele2 == bim_a2, "ok",
    Allele1 == bim_a2 & Allele2 == bim_a1, "flip",
    default = "drop")]
stats[allele_match == "flip", `:=`(Effect = -Effect, Freq1 = 1 - Freq1,
                                    Allele1 = bim_a1, Allele2 = bim_a2)]
n_before <- nrow(stats)
stats <- stats[!is.na(rsid) & allele_match != "drop"]
stats <- unique(stats, by = c("chr", "pos"))
message(sprintf("  %d UKB PDCD1 SNPs matched to bim (%d dropped)",
                nrow(stats), n_before - nrow(stats)))
## After this: stats$Allele1 == bim_a1 (same convention as T1D z-scores).

## Recompute loci on matched set
stats[, prev_pos  := shift(pos, 1L, type = "lag"), by = chr]
stats[, new_locus := is.na(prev_pos) | (pos - prev_pos) > gap_mb * 1e6]
stats[, locus_n   := cumsum(new_locus), by = chr]
stats[, locus_id  := paste0("chr", chr, "_locus", locus_n)]
stats[, c("prev_pos", "new_locus", "locus_n") := NULL]
message(sprintf("  %d loci after bim matching", uniqueN(stats$locus_id)))

## Remap top_snps_ukb locus_id to post-match locus_id
top_snps_ukb <- merge(top_snps_ukb,
                       stats[, .(chr, pos, locus_id_matched = locus_id)],
                       by.x = c("chr", "top_pos"), by.y = c("chr", "pos"),
                       all.x = FALSE)
top_snps_ukb[, locus_id := locus_id_matched][, locus_id_matched := NULL]

locus_n_alpha <- stats[, .(n_snps_alpha = .N), by = locus_id]
checkpoint("after_bim_match")

## ---- Step 5: Open 1000G BED and process each locus ----
## ---- Step 5: Build JSON with T1D z-scores and run zarr LD on genoscores ----
message("Building per-locus JSON (alpha + gamma) for genoscores...")
locus_ids <- unique(stats$locus_id)

job <- list(
    zarr_dir     = zarr_dir,
    min_eig_frac = min_eig_frac,
    N_gamma      = N_gamma,
    p_case       = p_case,
    loci = lapply(locus_ids, function(lid) {
        alpha_snps      <- stats[locus_id == lid]
        candidate_rsids <- intersect(alpha_snps$rsid, t1d$rsid)
        if (length(candidate_rsids) == 0L) return(NULL)
        alpha_snps <- alpha_snps[rsid %in% candidate_rsids]
        z_vec      <- setNames(t1d[.(candidate_rsids)]$z_t1d, candidate_rsids)
        list(
            locus_id = lid,
            chr      = alpha_snps$chr[1L],
            snps     = lapply(seq_len(nrow(alpha_snps)), function(i) list(
                rsid            = alpha_snps$rsid[i],
                Allele1         = alpha_snps$Allele1[i],
                Allele2         = alpha_snps$Allele2[i],
                Effect          = alpha_snps$Effect[i],
                Freq1           = alpha_snps$Freq1[i],
                StdErr          = alpha_snps$StdErr[i],
                TotalSampleSize = alpha_snps$TotalSampleSize[i],
                pos             = alpha_snps$pos[i],
                z_gamma         = z_vec[alpha_snps$rsid[i]]
            ))
        )
    })
)
job$loci <- Filter(Negate(is.null), job$loci)

json_in  <- tempfile(fileext = ".json")
json_err <- tempfile(fileext = ".err")
json_out <- tempfile(fileext = ".json")
write_json(job, json_in, auto_unbox = TRUE, digits = 8)
on.exit({ unlink(c(json_in, json_err, json_out)) }, add = TRUE)

message(sprintf("Sending %d loci to genoscores via SSH...", length(job$loci)))
ret <- system(sprintf(
    "scp -q run_zarr_loci.py '%s':~/ && ssh '%s' 'python3 ~/run_zarr_loci.py' < '%s' > '%s' 2>'%s'",
    genoscores_host, genoscores_host, json_in, json_out, json_err
))
if (file.exists(json_err) && file.size(json_err) > 0L)
    message(paste(readLines(json_err, warn = FALSE), collapse = "\n"))
if (ret != 0L) stop("run_zarr_loci.py failed on genoscores (exit ", ret, ")")
checkpoint("after_loci")

response <- read_json(json_out, simplifyVector = FALSE)$loci
gamma_dt <- rbindlist(lapply(response, function(r) {
    if (is.null(r$gamma_hat)) return(NULL)
    data.table(
        qtlname      = r$locus_id,
        chr          = as.integer(r$chr),
        locus_start  = as.integer(r$locus_start),
        locus_end    = as.integer(r$locus_end),
        n_snps_gamma = as.integer(r$n_snps),
        eff_rank     = as.integer(r$eff_rank),
        gamma_hat    = as.numeric(r$gamma_hat),
        se.gamma_hat = as.numeric(r$se_gamma_hat)
    )
}))
setorder(gamma_dt, chr, locus_start)

gamma_dt <- locus_n_alpha[gamma_dt, on = .(locus_id = qtlname)]
setnames(gamma_dt, "locus_id", "qtlname")
setcolorder(gamma_dt, c("qtlname", "chr", "locus_start", "locus_end",
                         "n_snps_alpha", "n_snps_gamma", "eff_rank",
                         "gamma_hat", "se.gamma_hat"))

## ---- Step 6: Report and save ----  (formerly Step 6)
message(sprintf("\nCompleted: %d loci with gamma estimates (of %d PDCD1 loci)",
                nrow(gamma_dt), length(locus_ids)))
print(gamma_dt[, .(qtlname, chr, n_snps_alpha, n_snps_gamma, eff_rank,
                    gamma_hat, se.gamma_hat)])

saveRDS(gamma_dt, "gamma_t1d_ukb.rds")
message("Saved: gamma_t1d_ukb.rds")

## ---- Step 7: Nearest protein-coding gene to top SNP per locus ----
message("Looking up nearest protein-coding gene per locus...")
genes <- fread("refGene.txt.gz", header = FALSE,
               select = c(2L, 3L, 5L, 6L, 13L),
               col.names = c("accession", "chrom", "txStart", "txEnd", "symbol"))
genes[, chr_g := suppressWarnings(as.integer(sub("^chr", "", chrom)))]
genes <- genes[!is.na(chr_g) & startsWith(accession, "NM_")]
genes <- genes[, .(txStart = min(txStart), txEnd = max(txEnd)),
                by = .(chr = chr_g, symbol)]

top_snps_ukb <- top_snps_ukb[locus_id %in% gamma_dt$qtlname]

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
