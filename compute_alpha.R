## compute_alpha.R
## Compute instrument-exposure coefficient alpha_hat for each locus
## from SCALLOP/UKBB meta-analysis summary statistics for PDCD1.
##
## LD source: 1000G EUR zarr store on genoscores (503 samples, GRCh37).
## Replaces the previous BEDMatrix approach; uses the same
## genoscores_zarr_loci.py pipeline as compute_alpha_ukb.R.

library(data.table)
library(jsonlite)
source("ld_functions.R")

## ---- LD reference panel paths (on genoscores) ----
zarr_dir_ukb <- "/opt/datastore/genome/LD_Eur_18mvariants/int8"
zarr_dir_1kg <- "/opt/datastore/genome/1000G/1kg_ld_eur/ld_1kg_eur"

## ---- Parameters ----
stats_file        <- "pdcd1_stats/OID00791_OID21396_SCALLOP_UKBB_MA_rsidannotated_filtered.tsv"
zarr_dir          <- zarr_dir_1kg   # 1000G EUR zarr
genoscores_host   <- "pmckeigue@genoscores.cphs.mvm.ed.ac.uk"
gap_mb            <- 1.0
window_mb         <- 1.0
min_eig_frac      <- 0.01
p_hit             <- 1e-6
p_cand            <- 1e-5
exclude_hla       <- TRUE
min_maf_threshold <- 0.02
out_file          <- "alpha_pdcd1_OID00791.rds"

## ---- Step 1: Load summary statistics ----
## SCALLOP/UKBB stats include rsid directly — no bim coordinate lookup needed.
message("Loading SCALLOP/UKBB meta-analysis summary statistics...")
stats <- fread(stats_file,
               select = c("chr","pos","rsid","Allele1","Allele2",
                          "Freq1","Effect","StdErr","TotalSampleSize","Pvalue"))
stats[, Allele1 := toupper(Allele1)]
stats[, Allele2 := toupper(Allele2)]
setorder(stats, chr, pos)
message(sprintf("  %d SNPs loaded", nrow(stats)))

is_ambiguous <- function(x, y)
    (x=="A"&y=="T")|(x=="T"&y=="A")|(x=="C"&y=="G")|(x=="G"&y=="C")
n_before <- nrow(stats)
stats <- stats[!is_ambiguous(Allele1, Allele2)]
message(sprintf("  %d SNPs after dropping strand-ambiguous (%d dropped)",
                nrow(stats), n_before - nrow(stats)))
stats <- stats[pmin(Freq1, 1 - Freq1) >= min_maf_threshold]
message(sprintf("  %d SNPs after MAF >= %.2f filter", nrow(stats), min_maf_threshold))

## ---- Step 2: Define loci ----
message("Defining loci (hit-anchored expansion)...")
stats <- define_loci_expanded(stats, p_hit=p_hit, p_cand=p_cand,
                               window_mb=window_mb, gap_mb=gap_mb,
                               exclude_hla=exclude_hla)
message(sprintf("  %d SNPs in %d loci", nrow(stats), uniqueN(stats$locus_id)))

top_snps <- stats[stats[, .I[which.min(Pvalue)], by=locus_id]$V1,
                  .(locus_id, chr, top_pos=pos)]
stats[, Pvalue := NULL]

## ---- Step 3: Build JSON and run zarr LD on genoscores ----
locus_ids       <- unique(stats$locus_id)
single_in_stats <- stats[, .N, by=locus_id][N == 1L, locus_id]

job <- list(
    zarr_dir     = zarr_dir,
    min_eig_frac = min_eig_frac,
    loci = lapply(locus_ids, function(lid) {
        snps <- stats[locus_id == lid]
        list(locus_id = lid, chr = snps$chr[1L],
             snps = lapply(seq_len(nrow(snps)), function(i) list(
                 rsid            = snps$rsid[i],
                 Allele1         = snps$Allele1[i],
                 Allele2         = snps$Allele2[i],
                 Effect          = snps$Effect[i],
                 Freq1           = snps$Freq1[i],
                 StdErr          = snps$StdErr[i],
                 TotalSampleSize = snps$TotalSampleSize[i],
                 pos             = snps$pos[i])))
    })
)

json_in  <- tempfile(fileext=".json")
json_err <- tempfile(fileext=".err")
json_out <- tempfile(fileext=".json")
write_json(job, json_in, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_in, json_err, json_out)), add=TRUE)

message(sprintf("Sending %d loci to genoscores via SSH (1000G EUR zarr)...",
                length(locus_ids)))
ret <- system(sprintf(
    "scp -q genoscores_zarr_loci.py '%s':~/ && ssh '%s' 'python3 ~/genoscores_zarr_loci.py' < '%s' > '%s' 2>'%s'",
    genoscores_host, genoscores_host, json_in, json_out, json_err))
if (file.exists(json_err) && file.size(json_err) > 0L)
    message(paste(readLines(json_err, warn=FALSE), collapse="\n"))
if (ret != 0L) stop("genoscores_zarr_loci.py failed (exit ", ret, ")")

response <- read_json(json_out, simplifyVector=FALSE)$loci
alpha_dt <- rbindlist(lapply(response, function(r) data.table(
    qtlname      = r$locus_id,
    chr          = as.integer(r$chr),
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

locus_min_maf <- stats[, .(min_maf = min(pmin(Freq1, 1 - Freq1))), by=locus_id]
alpha_dt <- merge(alpha_dt, locus_min_maf, by.x="qtlname", by.y="locus_id", all.x=TRUE)

## ---- Step 4: Report and save ----
message(sprintf("\nCompleted: %d loci", nrow(alpha_dt)))
print(alpha_dt[, .(qtlname, chr, n_snps, eff_rank, alpha_hat, se.alpha_hat)])

single_snp <- alpha_dt[n_snps == 1L & qtlname %in% single_in_stats]
if (nrow(single_snp) > 0L) {
    single_stats <- stats[locus_id %in% single_snp$qtlname,
                          .(qtlname=locus_id, expected=abs(Effect))]
    check <- merge(single_snp[, .(qtlname, alpha_hat)], single_stats, by="qtlname")
    message(sprintf("Single-SNP check: max |alpha_hat - |Effect|| / |Effect| = %.2e",
                    check[, max(abs(alpha_hat - expected) / expected)]))
}

saveRDS(alpha_dt, out_file)
message(sprintf("Saved: %s", out_file))

## ---- Step 5: Nearest protein-coding gene to top SNP per locus ----
message("Looking up nearest protein-coding gene per locus...")
genes <- fread("refGene.txt.gz", header=FALSE,
               select=c(2L,3L,5L,6L,13L),
               col.names=c("accession","chrom","txStart","txEnd","symbol"))
genes[, chr_g := suppressWarnings(as.integer(sub("^chr","",chrom)))]
genes <- genes[!is.na(chr_g) & startsWith(accession,"NM_")]
genes <- genes[, .(txStart=min(txStart), txEnd=max(txEnd)), by=.(chr=chr_g, symbol)]

top_snps <- top_snps[locus_id %in% alpha_dt$qtlname]
find_nearest <- function(chr_val, pos_val) {
    g <- genes[chr == chr_val]
    if (nrow(g) == 0L) return(NA_character_)
    dist <- pmax(0L, pmax(g$txStart - pos_val, pos_val - g$txEnd))
    g$symbol[which.min(dist)]
}
top_snps[, gene := mapply(find_nearest, chr, top_pos)]
saveRDS(top_snps[, .(qtlname=locus_id, gene)], "locus_genes.rds")
message("Saved: locus_genes.rds")
