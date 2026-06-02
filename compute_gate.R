## compute_gate.R
## Compute the GATE (Genetically Aided Trait Estimation) score coefficient
## from summary-level per-pQTL statistics, and compare with the directly
## computed coefficient from individual-level SDRNT1BIO/GS data on diabepi.
##
## GATE score for individual i over trans-pQTLs:
##   G_i  = sum_{k in trans} alpha_hat_k * Z_{ik}
##   G_i* = G_i / SD_ref(G)    [SD from UKBB reference panel]
##
## Summary-level formula (derived from score test; see gate_derivation.Rmd):
##   beta_GATE  = SD_G * sum_k(alpha_k * gamma_k * V_k) / sum_k(alpha_k^2 * V_k)
##   SE(beta)   = SD_G / sqrt(sum_k alpha_k^2 * V_k) = 1 / sqrt(N * p * (1-p))
##   V_k        = 1 / SE(gamma_k)^2  [Fisher information for locus k]
##   SD_G       = sqrt(sum_k alpha_k^2 * sd_Z_k^2)
##
## Direct computation: sends combined GATE weights for all trans-pQTL SNPs to
## diabepi as a single composite-instrument job; diabepi accumulates
##   GATE_i = sum_{k,j} G_{ij} * (alpha_hat_k * beta_{m,kj})
## standardises, and regresses T1D on the score.

library(data.table)
library(jsonlite)
source("ld_functions.R")

## ---- Parameters ----
tar_file          <- "pdcd1_stats/PDCD1_Q15116_OID21396_v1_Oncology.tar"
bim_dir           <- "refpop/bim_by_chr"
zarr_dir          <- "/opt/datastore/genome/LD_Eur_18mvariants/int8"
genoscores_host   <- "pmckeigue@genoscores.cphs.mvm.ed.ac.uk"
diabepi_host      <- "pmckeigue@diabepi.igmm.ed.ac.uk"
gap_mb            <- 1.0;  window_mb      <- 1.0
min_eig_frac      <- 0.01; p_hit          <- 1e-6
p_cand            <- 1e-5; exclude_hla    <- TRUE
info_threshold    <- 0.3;  min_maf_threshold <- 0.01
N_cases           <- 4922L; N_ctrls <- 7452L
N_gamma           <- N_cases + N_ctrls
p_case            <- N_cases / N_gamma
pdcd1_mid         <- 241858000L   # PDCD1 chr2 hg38 midpoint

## ---- Part A: Summary-level GATE ----
message("=== Part A: Summary-level GATE ===")

alpha_dt <- readRDS("alpha_pdcd1_UKB.rds")
gamma_dt <- readRDS("gamma_t1d_gs.rds")

## Use sd.Z from gamma_dt (consistent with se.gamma_hat; both from the same zarr run).
## Using alpha's sd.Z for SD_G while V_k derives from gamma's se.gamma introduces
## a mismatch that makes the SE too small.  With gamma's sd.Z, the theoretical
## SE = 1/sqrt(N*p*(1-p)) is recovered exactly.
dt <- merge(
    alpha_dt[, .(qtlname, chr, locus_start, locus_end, alpha_hat)],
    gamma_dt[, .(qtlname, gamma_hat, se.gamma_hat, sd.Z)],
    by = "qtlname"
)
message(sprintf("  %d loci with both alpha and gamma estimates", nrow(dt)))

## Exclude cis-pQTL (PDCD1, chr2)
cis_locus <- dt[chr == "2", qtlname[which.min(abs(locus_start - pdcd1_mid))]]
message(sprintf("  Cis-pQTL excluded from GATE sum: %s", cis_locus))
tr <- dt[qtlname != cis_locus]
message(sprintf("  %d trans-pQTLs contribute to GATE score", nrow(tr)))

V_k  <- 1 / tr$se.gamma_hat^2
SD_G <- sqrt(sum(tr$alpha_hat^2 * tr$sd.Z^2))

num       <- sum(tr$alpha_hat * tr$gamma_hat * V_k)
denom     <- sum(tr$alpha_hat^2 * V_k)
beta_gate <- SD_G * num / denom
se_gate   <- SD_G / sqrt(denom)
z_gate    <- beta_gate / se_gate
p_gate    <- 2 * pnorm(-abs(z_gate))

se_gate_theory <- 1 / sqrt(N_gamma * p_case * (1 - p_case))

message(sprintf("  SD_G (reference panel):            %.4f", SD_G))
message(sprintf("  beta_GATE (summary-level):         %.4f", beta_gate))
message(sprintf("  SE (empirical):                    %.6f", se_gate))
message(sprintf("  SE (theory: 1/sqrt(Np(1-p))):      %.6f", se_gate_theory))
message(sprintf("  z_GATE:                            %.4f", z_gate))
message(sprintf("  p-value:                           %.4g", p_gate))

## ---- Part B: Direct GATE (genoscores + diabepi) ----
message("\n=== Part B: Direct GATE ===")

## Load and filter PDCD1 summary stats (same pipeline as compute_alpha_ukb.R)
message("Loading PDCD1 summary statistics...")
members <- system2("tar", c("-tf", tar_file), stdout = TRUE)
members <- members[grepl("\\.gz$", members)]
read_one <- function(m) {
    fread(cmd = sprintf("tar -xOf '%s' '%s' | zcat", tar_file, m),
          select = c("CHROM","GENPOS","ALLELE0","ALLELE1","A1FREQ",
                     "INFO","N","BETA","SE","LOG10P"),
          colClasses = list(character = c("ALLELE0","ALLELE1")))
}
stats <- rbindlist(lapply(members, read_one))
stats <- stats[LOG10P > 5]
setnames(stats,
         c("CHROM","GENPOS","ALLELE0","ALLELE1","A1FREQ","N","BETA","SE","LOG10P"),
         c("chr","pos","Allele2","Allele1","Freq1","TotalSampleSize","Effect","StdErr","log10p"))
stats[, Allele1 := toupper(Allele1)]; stats[, Allele2 := toupper(Allele2)]
stats[, chr := as.integer(chr)]; setorder(stats, chr, pos)
stats <- stats[INFO >= info_threshold]; stats[, INFO := NULL]
is_ambig <- function(x, y)
    (x=="A"&y=="T")|(x=="T"&y=="A")|(x=="C"&y=="G")|(x=="G"&y=="C")
stats <- stats[!is_ambig(Allele1, Allele2)]
stats <- stats[pmin(Freq1, 1 - Freq1) >= min_maf_threshold]
stats[, Pvalue := 10^(-log10p)]; stats[, log10p := NULL]

stats <- define_loci_expanded(stats, p_hit=p_hit, p_cand=p_cand,
                               window_mb=window_mb, gap_mb=gap_mb,
                               exclude_hla=exclude_hla)
bim_chrs <- unique(stats$chr)
stats <- rbindlist(lapply(bim_chrs, function(ch) {
    f <- file.path(bim_dir, sprintf("chr_%s.rds", ch))
    if (!file.exists(f)) return(NULL)
    bim_ch <- readRDS(f)[, .(chr, pos, rsid)]
    bim_ch[stats[chr == ch], on = c("chr","pos")]
}))
stats <- stats[!is.na(rsid)]
stats <- unique(stats, by = c("chr","pos"))
stats[, prev_pos  := shift(pos, 1L, type="lag"), by=chr]
stats[, new_locus := is.na(prev_pos) | (pos - prev_pos) > gap_mb * 1e6]
stats[, locus_n   := cumsum(new_locus), by=chr]
stats[, locus_id  := paste0("chr", chr, "_locus", locus_n)]
stats[, c("prev_pos","new_locus","locus_n","Pvalue") := NULL]
message(sprintf("  %d SNPs in %d loci after filtering + bim-match",
                nrow(stats), uniqueN(stats$locus_id)))

## Send alpha job to genoscores to get beta_m weights per locus
message("Sending alpha job to genoscores for beta_m weights...")
locus_ids <- unique(stats$locus_id)
alpha_job <- list(
    zarr_dir     = zarr_dir,
    min_eig_frac = min_eig_frac,
    loci = lapply(locus_ids, function(lid) {
        s <- stats[locus_id == lid]
        list(locus_id = lid, chr = s$chr[1L],
             snps = lapply(seq_len(nrow(s)), function(i) list(
                 rsid            = s$rsid[i],
                 Allele1         = s$Allele1[i],
                 Allele2         = s$Allele2[i],
                 Effect          = s$Effect[i],
                 Freq1           = s$Freq1[i],
                 StdErr          = s$StdErr[i],
                 TotalSampleSize = s$TotalSampleSize[i],
                 pos             = s$pos[i])))
    })
)
json_a_in  <- tempfile(fileext=".json")
json_a_out <- tempfile(fileext=".json")
json_a_err <- tempfile(fileext=".err")
write_json(alpha_job, json_a_in, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_a_in, json_a_out, json_a_err)), add=TRUE)

ret <- system(sprintf(
    "scp -q genoscores_zarr_loci.py '%s':~/ && ssh '%s' 'python3 ~/genoscores_zarr_loci.py' < '%s' > '%s' 2>'%s'",
    genoscores_host, genoscores_host, json_a_in, json_a_out, json_a_err))
if (file.exists(json_a_err) && file.size(json_a_err) > 0L)
    message(paste(readLines(json_a_err, warn=FALSE), collapse="\n"))
if (ret != 0L) stop("genoscores zarr alpha failed (exit ", ret, ")")

alpha_resp <- read_json(json_a_out, simplifyVector=FALSE)$loci

## Extract alpha_hat and beta_m_weights per locus
alpha_info <- rbindlist(lapply(alpha_resp, function(r) data.table(
    locus_id  = r$locus_id,
    alpha_hat = as.numeric(r$alpha_hat)
)))
weights_by_locus <- setNames(
    lapply(alpha_resp, function(r) r$beta_m_weights),
    sapply(alpha_resp, `[[`, "locus_id"))
message(sprintf("  Received alpha + weights for %d loci", nrow(alpha_info)))

## Restrict to exactly the same trans-pQTL loci used in Part A (summary-level)
trans_loci_b <- alpha_info[locus_id %in% tr$qtlname, locus_id]
message(sprintf("  %d trans-pQTL loci for GATE (same set as summary-level)",
                length(trans_loci_b)))

## Build combined GATE SNP list: weight_{kj} = alpha_hat_k * beta_m_kj
gate_snps <- rbindlist(lapply(trans_loci_b, function(lid) {
    w_list    <- weights_by_locus[[lid]]
    if (is.null(w_list) || length(w_list) == 0L) return(NULL)
    a_k       <- alpha_info[locus_id == lid, alpha_hat]
    rbindlist(lapply(w_list, function(w) data.table(
        rsid    = w$rsid,
        Allele1 = w$Allele1,
        weight  = a_k * as.numeric(w$weight)
    )))
}))

## Check for duplicate rsids across loci (should not occur for trans-pQTLs)
n_dup <- sum(duplicated(gate_snps$rsid))
if (n_dup > 0L) {
    warning(sprintf("%d duplicate rsids across trans-pQTL loci — summing weights", n_dup))
    gate_snps <- gate_snps[, .(Allele1 = Allele1[1L], weight = sum(weight)), by=rsid]
}

message(sprintf("  GATE composite instrument: %d SNPs across %d loci",
                nrow(gate_snps), length(trans_loci_b)))

## Build Allele2 lookup from stats
allele2_map <- unique(stats[rsid %in% gate_snps$rsid, .(rsid, Allele2)])
gate_snps   <- merge(gate_snps, allele2_map, by="rsid", all.x=TRUE)

## Send to diabepi as a single "GATE" composite-instrument locus
diabepi_job <- list(loci = list(list(
    locus_id = "GATE",
    chr      = 0L,
    snps     = lapply(seq_len(nrow(gate_snps)), function(i) list(
        rsid    = gate_snps$rsid[i],
        Allele1 = gate_snps$Allele1[i],
        Allele2 = if (is.na(gate_snps$Allele2[i])) "N" else gate_snps$Allele2[i],
        weight  = gate_snps$weight[i]
    ))
)))

json_d_in  <- tempfile(fileext=".json")
json_d_out <- tempfile(fileext=".json")
json_d_err <- tempfile(fileext=".err")
write_json(diabepi_job, json_d_in, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_d_in, json_d_out, json_d_err)), add=TRUE)

message("Sending GATE composite instrument to diabepi...")
ret <- system(sprintf(
    "scp -q diabepi_gwas_loci.R '%s':~/ && scp -q '%s' '%s':~/snp_job.json && ssh '%s' 'Rscript ~/diabepi_gwas_loci.R ~/snp_job.json' > '%s' 2>'%s'",
    diabepi_host, json_d_in, diabepi_host, diabepi_host, json_d_out, json_d_err))
if (file.exists(json_d_err) && file.size(json_d_err) > 0L)
    message(paste(readLines(json_d_err, warn=FALSE), collapse="\n"))
if (ret != 0L) stop("diabepi GATE job failed (exit ", ret, ")")

d_resp    <- read_json(json_d_out, simplifyVector=FALSE)
gate_locus <- d_resp$loci[[1L]]

beta_gate_direct <- as.numeric(gate_locus$gamma_direct)
se_gate_direct   <- as.numeric(gate_locus$se_gamma_direct)
z_gate_direct    <- beta_gate_direct / se_gate_direct
p_gate_direct    <- 2 * pnorm(-abs(z_gate_direct))

message(sprintf("  beta_GATE (direct):   %.4f", beta_gate_direct))
message(sprintf("  SE (direct):          %.6f", se_gate_direct))
message(sprintf("  z_GATE (direct):      %.4f", z_gate_direct))
message(sprintf("  p-value (direct):     %.4g", p_gate_direct))

## ---- Comparison and save ----
message("\n=== Comparison: summary-level vs direct ===")
message(sprintf("  z_GATE summary:  %.4f", z_gate))
message(sprintf("  z_GATE direct:   %.4f", z_gate_direct))
message(sprintf("  Difference:      %.4f", z_gate - z_gate_direct))
message(sprintf("  beta ratio (summary/direct, ~ SD_G_ref/SD_G_target): %.4f",
                beta_gate / beta_gate_direct))

gate_results <- data.table(
    method           = c("summary", "direct"),
    beta_GATE        = c(beta_gate, beta_gate_direct),
    se_GATE          = c(se_gate,   se_gate_direct),
    z_GATE           = c(z_gate,    z_gate_direct),
    p_value          = c(p_gate,    p_gate_direct),
    n_trans_pqtl     = c(nrow(tr),  length(trans_loci_b)),
    SD_G             = c(SD_G,      NA_real_)
)
print(gate_results)

saveRDS(gate_results, "gate_results.rds")
message("Saved: gate_results.rds")
