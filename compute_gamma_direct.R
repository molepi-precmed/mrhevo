## compute_gamma_direct.R
## Compare two methods for computing instrument-outcome coefficients gamma_hat:
##
##   (A) Analytical (zarr LD): per-SNP SNP-outcome coefficients + UKBB zarr LD
##       → gamma_hat via the scalar instrument projection formula
##
##   (B) Direct individual-level regression: compute scalar instrument score
##       Z_i = sum_j G_{ij} * beta_m_j for each individual, normalise to SD=1,
##       regress T1D case-control status on Z_i with covariates.
##
## Pipeline:
##   Step 1  Local: define loci, bim-match (same as compute_alpha_ukb.R)
##   Step 2  genoscores: alpha computation → alpha_hat, sd_Z, beta_m weights
##   Step 3  diabepi: per-SNP z-statistics + direct regression on Z_i
##   Step 4  genoscores: analytical gamma from per-SNP z-statistics
##   Step 5  Local: merge and save comparison

library(data.table)
library(jsonlite)
source("ld_functions.R")

## ---- Parameters ----
tar_file        <- "pdcd1_stats/PDCD1_Q15116_OID21396_v1_Oncology.tar"
z_scores_file   <- "z_scores_gs.rds"
bim_dir         <- "refpop/bim_by_chr"
zarr_dir        <- "/opt/datastore/genome/LD_Eur_18mvariants/int8"
genoscores_host <- "pmckeigue@genoscores.cphs.mvm.ed.ac.uk"
diabepi_host    <- "pmckeigue@diabepi.igmm.ed.ac.uk"
gap_mb          <- 1.0;  window_mb     <- 1.0
min_eig_frac    <- 0.01; p_hit         <- 1e-6
p_cand          <- 1e-5; exclude_hla   <- TRUE
info_threshold  <- 0.3
N_cases         <- 4922L; N_ctrls <- 7452L
N_gamma         <- N_cases + N_ctrls; p_case <- N_cases / N_gamma

.t0 <- proc.time()["elapsed"]

## ---- Step 1: Load PDCD1 stats, define loci, bim-match ----
message("Loading PDCD1 summary stats and bim-matching...")
members <- system2("tar", c("-tf", tar_file), stdout = TRUE)
members <- members[grepl("\\.gz$", members)]
read_one <- function(m) {
    dt <- fread(cmd = sprintf("tar -xOf '%s' '%s' | zcat", tar_file, m),
                select = c("CHROM","GENPOS","ALLELE0","ALLELE1","A1FREQ",
                           "INFO","N","BETA","SE","LOG10P"),
                colClasses = list(character = c("ALLELE0","ALLELE1")))
    dt[LOG10P > 5]
}
stats <- rbindlist(lapply(members, read_one))
setnames(stats,
         c("CHROM","GENPOS","ALLELE0","ALLELE1","A1FREQ","N","BETA","SE","LOG10P"),
         c("chr","pos","Allele2","Allele1","Freq1","TotalSampleSize","Effect","StdErr","log10p"))
stats[, Allele1 := toupper(Allele1)]; stats[, Allele2 := toupper(Allele2)]
stats[, chr := as.integer(chr)]; setorder(stats, chr, pos)
stats <- stats[INFO >= info_threshold]; stats[, INFO := NULL]
stats[, Pvalue := 10^(-log10p)]; stats[, log10p := NULL]
is_ambig <- function(x,y) (x=="A"&y=="T")|(x=="T"&y=="A")|(x=="C"&y=="G")|(x=="G"&y=="C")
stats <- stats[!is_ambig(Allele1, Allele2)]
stats <- define_loci_expanded(stats, p_hit=p_hit, p_cand=p_cand,
                               window_mb=window_mb, gap_mb=gap_mb,
                               exclude_hla=exclude_hla)
bim_chrs <- unique(stats$chr)
stats <- rbindlist(lapply(bim_chrs, function(ch) {
    f <- file.path(bim_dir, sprintf("chr_%s.rds", ch))
    if (!file.exists(f)) return(NULL)
    bim_ch <- readRDS(f)[, .(chr, pos, rsid, bim_a1=a1, bim_a2=a2)]
    bim_ch[stats[chr==ch], on=c("chr","pos")]
}))
stats[, allele_match := fcase(
    Allele1==bim_a1 & Allele2==bim_a2, "ok",
    Allele1==bim_a2 & Allele2==bim_a1, "flip", default="drop")]
stats[allele_match=="flip", `:=`(Effect=-Effect, Freq1=1-Freq1,
                                  Allele1=bim_a1, Allele2=bim_a2)]
stats <- stats[!is.na(rsid) & allele_match!="drop"]
stats <- unique(stats, by=c("chr","pos"))
stats[, prev_pos := shift(pos,1L,type="lag"), by=chr]
stats[, new_locus := is.na(prev_pos)|(pos-prev_pos)>gap_mb*1e6]
stats[, locus_n   := cumsum(new_locus), by=chr]
stats[, locus_id  := paste0("chr",chr,"_locus",locus_n)]
stats[, c("prev_pos","new_locus","locus_n") := NULL]
message(sprintf("  %d SNPs in %d loci", nrow(stats), uniqueN(stats$locus_id)))

## ---- Step 2: zarr alpha + beta_m weights (genoscores) ----
message("Computing alpha and beta_m weights via zarr LD on genoscores...")
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
json_alpha_in  <- tempfile(fileext=".json")
json_alpha_out <- tempfile(fileext=".json")
json_alpha_err <- tempfile(fileext=".err")
write_json(alpha_job, json_alpha_in, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_alpha_in, json_alpha_out, json_alpha_err)), add=TRUE)

ret <- system(sprintf(
    "scp -q run_zarr_loci.py '%s':~/ && ssh '%s' 'python3 ~/run_zarr_loci.py' < '%s' > '%s' 2>'%s'",
    genoscores_host, genoscores_host, json_alpha_in, json_alpha_out, json_alpha_err))
if (file.exists(json_alpha_err) && file.size(json_alpha_err)>0L)
    message(paste(readLines(json_alpha_err, warn=FALSE), collapse="\n"))
if (ret != 0L) stop("zarr alpha failed (exit ", ret, ")")
message(sprintf("  zarr alpha done: %.1fs", proc.time()["elapsed"] - .t0))

alpha_response <- read_json(json_alpha_out, simplifyVector=FALSE)$loci
## Extract alpha_hat, sd_Z and beta_m_weights per locus
alpha_info <- rbindlist(lapply(alpha_response, function(r) data.table(
    locus_id    = r$locus_id,
    alpha_hat   = as.numeric(r$alpha_hat),
    sd_Z        = as.numeric(r$sd_Z)
)))
weights_by_locus <- setNames(
    lapply(alpha_response, function(r) r$beta_m_weights),
    sapply(alpha_response, `[[`, "locus_id"))

## ---- Step 3: diabepi — per-SNP z-statistics + direct gamma ----
message("Running diabepi (z-statistics + direct composite regression)...")
## Intersect with available z-scores (same SNPs as before)
z_gs <- readRDS(z_scores_file)   # rsid, locus_id, z_GS (from previous run)
setkey(z_gs, rsid)

loci_with_weights <- names(weights_by_locus)
diabepi_job <- list(loci = lapply(loci_with_weights, function(lid) {
    s      <- stats[locus_id == lid]
    w_list <- weights_by_locus[[lid]]
    ## build rsid → weight lookup
    w_map  <- setNames(sapply(w_list, `[[`, "weight"),
                       sapply(w_list, `[[`, "rsid"))
    ## keep only SNPs that are both bim-matched AND have a weight
    s <- s[rsid %in% names(w_map)]
    if (nrow(s) == 0L) return(NULL)
    list(locus_id = lid, chr = s$chr[1L],
         snps = lapply(seq_len(nrow(s)), function(i) list(
             rsid    = s$rsid[i],
             Allele1 = s$Allele1[i],
             Allele2 = s$Allele2[i],
             weight  = w_map[s$rsid[i]])))
}))
diabepi_job$loci <- Filter(Negate(is.null), diabepi_job$loci)

json_d_in  <- tempfile(fileext=".json")
json_d_out <- tempfile(fileext=".json")
json_d_err <- tempfile(fileext=".err")
write_json(diabepi_job, json_d_in, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_d_in, json_d_out, json_d_err)), add=TRUE)

ret <- system(sprintf(
    "scp -q run_gwas_loci.R '%s':~/ && scp -q '%s' '%s':~/snp_job.json && ssh '%s' 'Rscript ~/run_gwas_loci.R ~/snp_job.json' > '%s' 2>'%s'",
    diabepi_host, json_d_in, diabepi_host, diabepi_host, json_d_out, json_d_err))
if (file.exists(json_d_err) && file.size(json_d_err)>0L)
    message(paste(readLines(json_d_err, warn=FALSE), collapse="\n"))
if (ret != 0L) stop("diabepi run failed (exit ", ret, ")")
message(sprintf("  diabepi done: %.1fs", proc.time()["elapsed"] - .t0))

d_response <- read_json(json_d_out, simplifyVector=FALSE)

## Extract per-SNP z-statistics and direct gamma
z_new_dt <- rbindlist(lapply(d_response$loci, function(l) {
    if (length(l$snps)==0L) return(NULL)
    data.table(locus_id = l$locus_id,
               rsid     = vapply(l$snps, `[[`, character(1), "rsid"),
               z_gamma  = vapply(l$snps, `[[`, double(1),    "z_gamma"))
}))

direct_dt <- rbindlist(lapply(d_response$loci, function(l) {
    if (is.null(l$gamma_direct)) return(NULL)
    data.table(locus_id       = l$locus_id,
               gamma_direct   = as.numeric(l$gamma_direct),
               se_gamma_direct = as.numeric(l$se_gamma_direct))
}))
message(sprintf("  direct gamma: %d loci", nrow(direct_dt)))

## ---- Step 4: zarr analytical gamma (genoscores) ----
message("Computing analytical gamma via zarr LD on genoscores...")
stats_z <- merge(stats, z_new_dt, by=c("locus_id","rsid"), all.x=FALSE)
locus_ids_z <- unique(stats_z$locus_id)

gamma_job <- list(
    zarr_dir     = zarr_dir,
    min_eig_frac = min_eig_frac,
    N_gamma      = N_gamma,
    p_case       = p_case,
    loci = lapply(locus_ids_z, function(lid) {
        s <- stats_z[locus_id == lid]
        list(locus_id = lid, chr = s$chr[1L],
             snps = lapply(seq_len(nrow(s)), function(i) list(
                 rsid            = s$rsid[i],
                 Allele1         = s$Allele1[i],
                 Allele2         = s$Allele2[i],
                 Effect          = s$Effect[i],
                 Freq1           = s$Freq1[i],
                 StdErr          = s$StdErr[i],
                 TotalSampleSize = s$TotalSampleSize[i],
                 pos             = s$pos[i],
                 z_gamma         = s$z_gamma[i])))
    })
)
json_g_in  <- tempfile(fileext=".json")
json_g_out <- tempfile(fileext=".json")
json_g_err <- tempfile(fileext=".err")
write_json(gamma_job, json_g_in, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_g_in, json_g_out, json_g_err)), add=TRUE)

ret <- system(sprintf(
    "ssh '%s' 'python3 ~/run_zarr_loci.py' < '%s' > '%s' 2>'%s'",
    genoscores_host, json_g_in, json_g_out, json_g_err))
if (file.exists(json_g_err) && file.size(json_g_err)>0L)
    message(paste(readLines(json_g_err, warn=FALSE), collapse="\n"))
if (ret != 0L) stop("zarr gamma failed (exit ", ret, ")")
message(sprintf("  zarr gamma done: %.1fs", proc.time()["elapsed"] - .t0))

g_response <- read_json(json_g_out, simplifyVector=FALSE)$loci
analytic_dt <- rbindlist(lapply(g_response, function(r) {
    if (is.null(r$gamma_hat)) return(NULL)
    data.table(locus_id     = r$locus_id,
               gamma_analytic    = as.numeric(r$gamma_hat),
               se_gamma_analytic = as.numeric(r$se_gamma_hat),
               sd_Z         = as.numeric(r$sd_Z))
}))

## ---- Step 5: Merge and save ----
comp <- merge(analytic_dt, direct_dt, by="locus_id")

## Load nearest gene for labelling
locus_genes <- readRDS("locus_genes_ukb.rds")
comp <- merge(comp, locus_genes, by.x="locus_id", by.y="qtlname", all.x=TRUE)

message(sprintf("\nCompleted: %d loci with both analytical and direct gamma",
                nrow(comp)))

## Scale both to per-SD(Z) units
comp[, gamma_analytic_sdz := gamma_analytic * sd_Z]
comp[, gamma_direct_sdz   := gamma_direct   * sd_Z]

r_cor <- cor(comp$gamma_analytic_sdz, comp$gamma_direct_sdz, use="complete.obs")
message(sprintf("Pearson r (analytical vs direct, per-SD(Z) units): %.4f", r_cor))
print(comp[, .(locus_id, gene, gamma_analytic_sdz, gamma_direct_sdz,
               se_gamma_analytic, se_gamma_direct)])

saveRDS(comp, "gamma_comparison_analytic_vs_direct.rds")
message("Saved: gamma_comparison_analytic_vs_direct.rds")
message(sprintf("Total elapsed: %.1fs", proc.time()["elapsed"] - .t0))
