#!/usr/bin/env Rscript
## diabepi_gwas_loci.R
## Runs on diabepi; master copy kept locally under git.
## Reads JSON from stdin, writes JSON to stdout.
##
## Input JSON:
##   { "loci": [{ "locus_id": "chr2_locus1", "chr": 2,
##               "snps": [{ "rsid": "rs...", "Allele1": "A", "Allele2": "G" }] }] }
##
## Output JSON:
##   { "N_cases": 4922, "N_ctrls": 7452,
##     "loci": [{ "locus_id": "chr2_locus1",
##               "snps": [{ "rsid": "rs...", "z_gamma": 1.23 }] }] }
##
## z_gamma is the Wald z-statistic from logistic regression of T1D on genotype,
## adjusted for sex, age2016, PC1, PC2, PC3; signed in the Allele1 direction.
## No individual-level data is written anywhere.

library(data.table)
library(jsonlite)
library(BEDMatrix)

bed_prefix <- "/opt/shared/project/type1bio/data/2014/2016-11-24_V22/plink/GS_T1biochrall_filtered_rsmapped"
pheno_file <- paste0(
    "/opt/shared/project/type1bio/data/gwas_analysis/gen_scotland/phenotypes/",
    "PHENO_T1DM_GST1B_AS_NOREL_NOMODY_NOPT2_NOCPEPAB_NONEO_NORELPCA.samples")

## ---- Open plink files ----
message("Opening BED file (memory-mapped)...")
bed <- BEDMatrix(bed_prefix, simple_names = TRUE)   # rows = individuals, cols = SNPs
bim <- fread(paste0(bed_prefix, ".bim"), header = FALSE,
             col.names = c("chr", "rsid", "cm", "pos", "A1", "A2"))
fam <- fread(paste0(bed_prefix, ".fam"), header = FALSE,
             col.names = c("FID", "IID", "PID", "MID", "sex_fam", "pheno_fam"))
message(sprintf("  BED: %d individuals x %d SNPs", nrow(fam), nrow(bim)))

## ---- Read phenotype file ----
## .samples format: line 1 = col names, line 2 = type codes, lines 3+ = data
pheno_raw <- fread(pheno_file, skip = 2L, header = FALSE,
                   col.names = c("ID_1", "ID_2", "missing",
                                 "PHENO", "sex", "age2016", "PC1", "PC2", "PC3"))
pheno_raw[, PHENO    := as.integer(PHENO)]
pheno_raw[, sex      := as.numeric(sex)]
pheno_raw[, age2016  := as.numeric(age2016)]
message(sprintf("  Phenotype file: %d records", nrow(pheno_raw)))

## ---- Sample intersection ----
m <- match(fam$IID, pheno_raw$ID_1)
n_matched <- sum(!is.na(m))
message(sprintf("  FAM IID matched to phenotype: %d / %d (%.0f%%)",
                n_matched, nrow(fam), 100 * n_matched / nrow(fam)))
if (n_matched / nrow(fam) < 0.5)
    stop(sprintf(
        "Only %d/%d FAM IIDs matched phenotype IDs — ID scheme mismatch, check plink files",
        n_matched, nrow(fam)))
pheno_matched <- pheno_raw[m, ]   # NA rows where no phenotype record

## ---- Complete cases ----
complete <- !is.na(pheno_matched$PHENO) &
            pheno_matched$PHENO %in% 0:1 &
            !is.na(pheno_matched$sex)    &
            !is.na(pheno_matched$age2016) &
            !is.na(pheno_matched$PC1)
sample_idx <- which(complete)
pheno_vec  <- pheno_matched$PHENO[complete]
cov_df     <- pheno_matched[complete, .(sex, age2016, PC1, PC2, PC3)]

N_cases <- sum(pheno_vec == 1L)
N_ctrls <- sum(pheno_vec == 0L)
message(sprintf("  Complete cases: %d cases, %d controls", N_cases, N_ctrls))

## ---- Build BIM rsid index ----
bim_col_idx <- setNames(seq_len(nrow(bim)), bim$rsid)

## ---- Fit null logistic model (covariates only) ----
message("Fitting null logistic model (covariates only)...")
null_df  <- data.frame(pheno = pheno_vec, cov_df)
null_fit <- glm(pheno ~ ., data = null_df, family = binomial)
mu_null  <- fitted(null_fit)
r_null   <- pheno_vec - mu_null     # score residuals under null
w_null   <- mu_null * (1 - mu_null) # Fisher information weights

## Precompute Cholesky of Z'WZ for partial-information correction.
## Without this, V = x'Wx overestimates information because it includes
## the contribution of x's mean (captured by the intercept) and any
## correlation of x with covariates; the corrected V removes these.
Z_mat     <- model.matrix(null_fit)        # n × 6 (intercept + 5 covariates)
WZ_mat    <- Z_mat * w_null                # n × 6, rows scaled by w_null
ZtWZ      <- crossprod(WZ_mat, Z_mat)      # 6 × 6
chol_ZtWZ <- chol(ZtWZ)                   # upper Cholesky: t(U) %*% U = ZtWZ
rm(null_df, null_fit, Z_mat)
message("  Null model fitted.")

## ---- Process loci ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript diabepi_gwas_loci.R <input.json>")
job <- fromJSON(args[1L], simplifyVector = FALSE)
message(sprintf("Processing %d loci...", length(job$loci)))

locus_results <- lapply(job$loci, function(locus) {
    lid      <- locus$locus_id
    snps_in  <- locus$snps
    rsids    <- vapply(snps_in, `[[`, character(1), "rsid")
    allele1  <- setNames(vapply(snps_in, `[[`, character(1), "Allele1"), rsids)

    col_j <- bim_col_idx[rsids]           # NA for rsids not in BIM
    n_found <- sum(!is.na(col_j))
    if (n_found == 0L) {
        message(sprintf("  %s: no rsids found in BIM, skipping", lid))
        return(list(locus_id = lid, chr = locus$chr, snps = list()))
    }
    message(sprintf("  %s: %d/%d rsids found in BIM", lid, n_found, length(rsids)))

    snp_out <- lapply(seq_along(rsids), function(k) {
        rsid <- rsids[k]
        j    <- col_j[k]
        if (is.na(j)) return(NULL)

        geno <- tryCatch(as.integer(bed[sample_idx, j]), error = function(e) NULL)
        if (is.null(geno)) return(NULL)

        ## Allele alignment: BEDMatrix returns count of BIM A1
        ## If BIM A1 != Allele1, flip to Allele1 direction
        bim_a1 <- toupper(bim$A1[j])
        req_a1 <- toupper(allele1[rsid])
        if (bim_a1 != req_a1)
            geno <- 2L - geno

        keep <- !is.na(geno)
        if (sum(keep) < 100L) return(NULL)   # too many missing

        x <- as.double(geno[keep])
        U <- sum(x * r_null[keep])

        ## Partial information: V_raw - projection onto null-model covariate space.
        ## Corrects for x's correlation with the intercept (allele frequency mean)
        ## and with sex, age, PCs.
        V_raw  <- sum(x * x * w_null[keep])
        ZtWx   <- crossprod(WZ_mat[keep, , drop = FALSE], x)  # 6-vector
        c_vec  <- forwardsolve(t(chol_ZtWZ), ZtWx)
        V      <- V_raw - sum(c_vec * c_vec)
        if (V <= 0) return(NULL)

        list(rsid     = rsid,
             z_gamma  = U / sqrt(V),
             gamma    = U / V,
             se_gamma = 1 / sqrt(V))
    })
    snp_out <- Filter(Negate(is.null), snp_out)
    result  <- list(locus_id = lid, chr = locus$chr, snps = snp_out)

    ## ---- Direct regression of T1D on composite instrument score ----
    ## Only run when instrument weights (beta_m in Allele1 direction) are provided.
    has_weights <- !is.null(snps_in[[1L]]$weight)
    if (has_weights) {
        weights <- setNames(vapply(snps_in, `[[`, double(1), "weight"), rsids)

        ## Accumulate Z_i = sum_j geno_j * weight_j with mean imputation for missing
        Z_score <- rep(0, length(sample_idx))
        n_contrib <- 0L
        for (k in seq_along(rsids)) {
            rsid <- rsids[k]; j <- col_j[k]
            if (is.na(j)) next
            geno <- tryCatch(as.integer(bed[sample_idx, j]), error = function(e) NULL)
            if (is.null(geno)) next
            bim_a1 <- toupper(bim$A1[j])
            req_a1 <- toupper(allele1[rsid])
            if (bim_a1 != req_a1) geno <- 2L - geno
            ## mean-impute missing genotypes before accumulating
            mu <- mean(geno, na.rm = TRUE)
            geno[is.na(geno)] <- mu
            Z_score   <- Z_score + geno * weights[rsid]
            n_contrib <- n_contrib + 1L
        }

        if (n_contrib > 0L && sd(Z_score) > 0) {
            Z_scaled <- Z_score / sd(Z_score)   # unit-SD instrument in target dataset

            ## Score test (same partial-information correction as per-SNP code)
            U_z    <- sum(Z_scaled * r_null)
            V_raw_z <- sum(Z_scaled * Z_scaled * w_null)
            ZtWz   <- crossprod(WZ_mat, Z_scaled)
            c_vec_z <- forwardsolve(t(chol_ZtWZ), ZtWz)
            V_z    <- V_raw_z - sum(c_vec_z * c_vec_z)
            if (V_z > 0) {
                result$gamma_direct    <- U_z / V_z
                result$se_gamma_direct <- 1 / sqrt(V_z)
                message(sprintf("  %s: gamma_direct = %.4f (SE %.4f) [score test]",
                                lid, result$gamma_direct, result$se_gamma_direct))
            }

        }
    }
    result
})

## ---- Write output ----
out <- list(N_cases = N_cases, N_ctrls = N_ctrls, loci = locus_results)
write_json(out, stdout(), auto_unbox = TRUE, digits = 8)
message("Done.")
