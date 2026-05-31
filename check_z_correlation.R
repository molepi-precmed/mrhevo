## check_z_correlation.R
## Compare per-SNP z-scores from Iakovliev et al. (gwasresults.RData.gz)
## with z-scores computed from individual-level SDRNT1BIO/GS genotypes.
## Checks allele coding alignment and reports Pearson correlation.
## Saves z_scores_gs.rds (per-SNP GS z-scores) for future use.

library(data.table)
library(jsonlite)
library(ggplot2)
library(ggrepel)
source("ld_functions.R")

## ---- Parameters (same as compute_gamma_indiv.R) ----
tar_file       <- "pdcd1_stats/PDCD1_Q15116_OID21396_v1_Oncology.tar"
t1d_file       <- "t1d_stats/data/gwasresults.RData.gz"
bim_dir        <- "refpop/bim_by_chr"
diabepi_host   <- "pmckeigue@diabepi.igmm.ed.ac.uk"
gap_mb         <- 1.0; window_mb <- 1.0
min_eig_frac   <- 0.01; p_hit <- 1e-6; p_cand <- 1e-5
exclude_hla    <- TRUE; info_threshold <- 0.3

## ---- Step 1: Load and filter PDCD1 summary stats ----
message("Loading PDCD1 summary stats...")
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

## ---- Step 2: Define loci and bim-match ----
message("Defining loci and bim-matching...")
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
stats[, locus_n := cumsum(new_locus), by=chr]
stats[, locus_id := paste0("chr",chr,"_locus",locus_n)]
stats[, c("prev_pos","new_locus","locus_n") := NULL]
message(sprintf("  %d SNPs in %d loci after bim-matching", nrow(stats), uniqueN(stats$locus_id)))

## ---- Step 3: Load Iakovliev z-scores ----
message("Loading Iakovliev T1D z-scores...")
tmp_rdata <- tempfile(fileext=".RData")
system2("gunzip", args=c("-c", t1d_file), stdout=tmp_rdata)
load(tmp_rdata); unlink(tmp_rdata)
t1d <- as.data.table(results); rm(results)
t1d[, CHR := as.integer(as.character(CHR))]
setnames(t1d, c("SNP","CHR","z"), c("rsid","chr","z_IA"))
t1d[, c("BP","minuslog10p") := NULL]
if (anyDuplicated(t1d$rsid))
    t1d <- t1d[t1d[, .I[which.max(abs(z_IA))], by=rsid]$V1]
setkey(t1d, rsid)
message(sprintf("  %d Iakovliev SNPs loaded", nrow(t1d)))

## Restrict PDCD1 locus SNPs to those present in Iakovliev
stats_ia <- stats[rsid %in% t1d$rsid]
message(sprintf("  %d locus SNPs present in Iakovliev data", nrow(stats_ia)))
locus_ids <- unique(stats_ia$locus_id)

## ---- Step 4: Build JSON and run diabepi for GS z-scores ----
message(sprintf("Sending %d loci (%d SNPs) to diabepi...", length(locus_ids), nrow(stats_ia)))
snp_job <- list(loci = lapply(locus_ids, function(lid) {
    s <- stats_ia[locus_id==lid]
    list(locus_id=lid, chr=s$chr[1L],
         snps=lapply(seq_len(nrow(s)), function(i)
             list(rsid=s$rsid[i], Allele1=s$Allele1[i], Allele2=s$Allele2[i])))
}))
json_snp  <- tempfile(fileext=".json")
json_z    <- tempfile(fileext=".json")
json_zerr <- tempfile(fileext=".err")
write_json(snp_job, json_snp, auto_unbox=TRUE, digits=8)
on.exit(unlink(c(json_snp, json_z, json_zerr)), add=TRUE)

ret <- system(sprintf(
    "scp -q run_gwas_loci.R '%s':~/ && scp -q '%s' '%s':~/snp_job.json && ssh '%s' 'Rscript ~/run_gwas_loci.R ~/snp_job.json' > '%s' 2>'%s'",
    diabepi_host, json_snp, diabepi_host, diabepi_host, json_z, json_zerr))
if (file.exists(json_zerr) && file.size(json_zerr)>0L)
    message(paste(readLines(json_zerr, warn=FALSE), collapse="\n"))
if (ret!=0L) stop("run_gwas_loci.R failed on diabepi (exit ", ret, ")")

z_response <- read_json(json_z, simplifyVector=FALSE)
z_dt <- rbindlist(lapply(z_response$loci, function(l) {
    if (length(l$snps)==0L) return(NULL)
    data.table(locus_id = l$locus_id,
               rsid     = vapply(l$snps, `[[`, character(1), "rsid"),
               z_GS     = vapply(l$snps, `[[`, double(1),    "z_gamma"))
}))
message(sprintf("  Received z_GS for %d SNPs across %d loci", nrow(z_dt), uniqueN(z_dt$locus_id)))
saveRDS(z_dt, "z_scores_gs.rds")
message("Saved: z_scores_gs.rds")

## ---- Step 5: Merge and check allele coding ----
## Both z_IA and z_GS should be in the Allele1 (bim_a1) direction.
## z_GS: explicitly aligned to Allele1 by run_gwas_loci.R.
## z_IA: allele direction from Iakovliev; no flip was applied in compute_gamma_ukb.R,
##       so z_IA is in whatever direction Iakovliev used.
comp <- merge(z_dt, t1d[, .(rsid, z_IA)], by="rsid")
## Add Allele1 (bim_a1) for reference
comp <- merge(comp, stats[, .(rsid, locus_id, Allele1, bim_a1, bim_a2)], by=c("rsid","locus_id"))
message(sprintf("  %d SNPs matched in both datasets", nrow(comp)))

## Allele coding check: what fraction of z_IA and z_GS have the same sign?
same_sign_frac <- comp[, mean(sign(z_IA)==sign(z_GS))]
message(sprintf("  Fraction same sign (z_IA vs z_GS): %.3f", same_sign_frac))

if (same_sign_frac < 0.5) {
    message("  WARNING: Majority of z-scores have opposite signs — allele coding is flipped.")
    message("  Flipping z_IA for correlation computation.")
    comp[, z_IA := -z_IA]
} else {
    message("  Allele coding appears consistent (majority of z-scores have the same sign).")
}

## ---- Step 6: Pearson correlation ----
r <- comp[, cor(z_GS, z_IA, use="complete.obs")]
message(sprintf("\nPearson correlation (z_GS vs z_IA): r = %.4f", r))
message(sprintf("  N = %d SNPs", nrow(comp[!is.na(z_GS) & !is.na(z_IA)])))

## ---- Step 7: Scatter plot if r < 0.99 ----
if (abs(r) < 0.99) {
    message("  Correlation < 0.99 — saving scatter plot.")

    locus_label <- comp[, .(n=.N), by=locus_id][order(-n)]
    comp[, locus_short := sub("chr([^_]+)_locus.*", "chr\\1", locus_id)]

    p <- ggplot(comp, aes(x=z_GS, y=z_IA)) +
        geom_point(alpha=0.4, size=0.8) +
        geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") +
        annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.5,
                 label=sprintf("r = %.4f\nN = %d", r, nrow(comp)),
                 size=3.5) +
        labs(x = "z-score (GS individual-level logistic regression, bim_a1 direction)",
             y = "z-score (Iakovliev et al. GWAS summary statistics)",
             title = "Per-SNP T1D z-score comparison") +
        theme_bw()
    ggsave("z_score_comparison.pdf", p, width=7, height=6)
    message("Saved: z_score_comparison.pdf")
} else {
    message("  Correlation >= 0.99 — no scatter plot needed.")
}

message("\nDone.")
