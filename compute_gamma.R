## compute_gamma.R
## Compute instrument-outcome coefficient gamma_hat for each locus
## from GWAS summary statistics (z-scores) for SNP-T1D associations.
##
## LD source: 1000G EUR zarr store on genoscores (503 samples, GRCh37).
## Replaces the previous BEDMatrix approach; uses the same
## genoscores_zarr_loci.py pipeline as compute_gamma_ukb.R.
##
## Exposure: SCALLOP/UKBB meta-analysis summary stats for PDCD1 (OID00791).
## Outcome:  Iakovliev et al. 2023 T1D GWAS (AJHG, PMID 37164005).

library(data.table)
library(jsonlite)
source("ld_functions.R")

## ---- LD reference panel paths (on genoscores) ----
zarr_dir_ukb <- "/opt/datastore/genome/LD_Eur_18mvariants/int8"
zarr_dir_1kg <- "/opt/datastore/genome/1000G/1kg_ld_eur/ld_1kg_eur"

## ---- Parameters ----
stats_file        <- "pdcd1_stats/OID00791_OID21396_SCALLOP_UKBB_MA_rsidannotated_filtered.tsv"
t1d_file          <- "t1d_stats/data/gwasresults.RData.gz"
zarr_dir          <- zarr_dir_1kg   # 1000G EUR zarr
genoscores_host   <- "pmckeigue@genoscores.cphs.mvm.ed.ac.uk"
gap_mb            <- 1.0
window_mb         <- 1.0
min_eig_frac      <- 0.01
p_hit             <- 1e-6
p_cand            <- 1e-5
exclude_hla       <- TRUE
min_maf_threshold <- 0.02
out_file          <- "gamma_t1d.rds"

## T1D GWAS sample size (Iakovliev et al. 2023)
N_cases <- 4964L
N_ctrls <- 7497L
N_gamma <- N_cases + N_ctrls
p_case  <- N_cases / N_gamma

## ---- Step 1: Load T1D GWAS z-scores ----
message("Loading T1D GWAS summary statistics (Iakovliev et al.)...")
tmp_rdata <- tempfile(fileext=".RData")
system2("gunzip", args=c("-c", t1d_file), stdout=tmp_rdata)
load(tmp_rdata); unlink(tmp_rdata)
t1d <- as.data.table(results); rm(results)
t1d[, CHR := as.integer(as.character(CHR))]
setnames(t1d, c("SNP","CHR","z"), c("rsid","chr","z_t1d"))
t1d[, c("BP","minuslog10p") := NULL]
if (anyDuplicated(t1d$rsid))
    t1d <- t1d[t1d[, .I[which.max(abs(z_t1d))], by=rsid]$V1]
setkey(t1d, rsid)
message(sprintf("  %d T1D SNPs loaded", nrow(t1d)))

## ---- Step 2: Load PDCD1 summary statistics ----
## SCALLOP/UKBB stats include rsid directly — no bim coordinate lookup needed.
message("Loading SCALLOP/UKBB PDCD1 summary statistics...")
stats <- fread(stats_file,
               select=c("chr","pos","rsid","Allele1","Allele2",
                        "Freq1","Effect","StdErr","TotalSampleSize","Pvalue"))
stats[, Allele1 := toupper(Allele1)]
stats[, Allele2 := toupper(Allele2)]
setorder(stats, chr, pos)
message(sprintf("  %d PDCD1 SNPs loaded", nrow(stats)))

is_ambiguous <- function(x, y)
    (x=="A"&y=="T")|(x=="T"&y=="A")|(x=="C"&y=="G")|(x=="G"&y=="C")
stats <- stats[!is_ambiguous(Allele1, Allele2)]
stats <- stats[pmin(Freq1, 1 - Freq1) >= min_maf_threshold]
message(sprintf("  %d SNPs after strand-ambiguous and MAF >= %.2f filters",
                nrow(stats), min_maf_threshold))

## ---- Step 3: Define loci ----
message("Defining loci (hit-anchored expansion)...")
stats <- define_loci_expanded(stats, p_hit=p_hit, p_cand=p_cand,
                               window_mb=window_mb, gap_mb=gap_mb,
                               exclude_hla=exclude_hla)
message(sprintf("  %d PDCD1 SNPs in %d loci", nrow(stats), uniqueN(stats$locus_id)))
locus_n_alpha <- stats[, .(n_snps_alpha=.N), by=locus_id]
stats[, Pvalue := NULL]

## ---- Step 4: Build JSON (with T1D z-scores) and run zarr LD on genoscores ----
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
        list(locus_id = lid, chr = alpha_snps$chr[1L],
             snps = lapply(seq_len(nrow(alpha_snps)), function(i) list(
                 rsid            = alpha_snps$rsid[i],
                 Allele1         = alpha_snps$Allele1[i],
                 Allele2         = alpha_snps$Allele2[i],
                 Effect          = alpha_snps$Effect[i],
                 Freq1           = alpha_snps$Freq1[i],
                 StdErr          = alpha_snps$StdErr[i],
                 TotalSampleSize = alpha_snps$TotalSampleSize[i],
                 pos             = alpha_snps$pos[i],
                 z_gamma         = z_vec[alpha_snps$rsid[i]]
             )))
    })
)
job$loci <- Filter(Negate(is.null), job$loci)

json_in  <- tempfile(fileext=".json")
json_err <- tempfile(fileext=".err")
json_out <- tempfile(fileext=".json")
write_json(job, json_in, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_in, json_err, json_out)), add=TRUE)

message(sprintf("Sending %d loci to genoscores via SSH (1000G EUR zarr)...",
                length(job$loci)))
ret <- system(sprintf(
    "scp -q genoscores_zarr_loci.py '%s':~/ && ssh '%s' 'python3 ~/genoscores_zarr_loci.py' < '%s' > '%s' 2>'%s'",
    genoscores_host, genoscores_host, json_in, json_out, json_err))
if (file.exists(json_err) && file.size(json_err) > 0L)
    message(paste(readLines(json_err, warn=FALSE), collapse="\n"))
if (ret != 0L) stop("genoscores_zarr_loci.py failed (exit ", ret, ")")

response <- read_json(json_out, simplifyVector=FALSE)$loci
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

gamma_dt <- locus_n_alpha[gamma_dt, on=.(locus_id=qtlname)]
setnames(gamma_dt, "locus_id", "qtlname")
setcolorder(gamma_dt, c("qtlname","chr","locus_start","locus_end",
                         "n_snps_alpha","n_snps_gamma","eff_rank",
                         "gamma_hat","se.gamma_hat"))

## ---- Step 5: Report and save ----
message(sprintf("\nCompleted: %d loci with gamma estimates (of %d PDCD1 loci)",
                nrow(gamma_dt), length(locus_ids)))
print(gamma_dt[, .(qtlname, chr, n_snps_alpha, n_snps_gamma, eff_rank,
                    gamma_hat, se.gamma_hat)])

saveRDS(gamma_dt, out_file)
message(sprintf("Saved: %s", out_file))
