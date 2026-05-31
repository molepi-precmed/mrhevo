## ld_functions.R
## LD matrix access and LD-adjusted scalar instrument construction.
##
## Supports two LD sources:
##   (a) UKBB EUR Zarr store via magenpy/Python (get_locus_ld)
##   (b) 1000G EUR BED file via BEDMatrix (get_locus_ld_1kg)
##
## Designed so this script can eventually be hosted on a dedicated LD server:
##   - Client prepares per-locus SNP tables and z-score vectors
##   - Server runs locus_compute_alpha() / locus_compute_gamma() with LD matrix access
##   - Server returns coefficient lists; client assembles output tables

library(data.table)

## ---- UKBB Zarr LD (Python/magenpy) ----

library(reticulate)
source_python("get_ld_matrix.py")

## ---- 1000G BED LD (BEDMatrix) ----

library(BEDMatrix)

## Per-chromosome bim cache (avoids re-reading per locus call).
.bim_cache <- new.env(hash = TRUE, parent = emptyenv())

.load_chr_bim <- function(bim_dir, chrom) {
    key <- as.character(chrom)
    if (!exists(key, envir = .bim_cache)) {
        f <- file.path(bim_dir, sprintf("chr_%s.rds", chrom))
        dt <- readRDS(f)
        setkey(dt, rsid)
        assign(key, dt, envir = .bim_cache)
    }
    get(key, envir = .bim_cache)
}

## Extract dense LD correlation submatrix from 1000G BED file.
##
## @param bed     BEDMatrix object (opened once by caller).
## @param bim_dir Path to directory of per-chromosome bim RDS files (with global_idx).
## @param chrom   Chromosome number (integer).
## @param rsids   Character vector of rsIDs.
##
## @return list(C_g, snp_df) — same format as get_locus_ld() — or list(NULL, NULL).
get_locus_ld_1kg <- function(bed, bim_dir, chrom, rsids) {
    f <- file.path(bim_dir, sprintf("chr_%s.rds", chrom))
    if (!file.exists(f)) return(list(NULL, NULL))
    chr_bim <- .load_chr_bim(bim_dir, chrom)
    matched <- chr_bim[rsid %in% rsids]
    if (nrow(matched) == 0L) return(list(NULL, NULL))

    G <- bed[, matched$global_idx, drop = FALSE]   # 503 × k
    G_std <- scale(G)

    ## Drop monomorphic or all-NA columns
    ok <- colSums(!is.na(G_std)) > 1L & !is.nan(G_std[1L, ])
    if (!any(ok)) return(list(NULL, NULL))
    G_std   <- G_std[, ok, drop = FALSE]
    matched <- matched[ok]

    C_g <- cor(G_std, use = "pairwise.complete.obs")
    ## Replace any residual NA (shouldn't occur with 503 samples) with 0
    C_g[is.na(C_g)] <- 0

    snp_df <- data.frame(SNP = matched$rsid,
                         A1  = matched$a1,
                         A2  = matched$a2,
                         POS = matched$pos,
                         stringsAsFactors = FALSE)
    list(C_g, snp_df)
}

## ---- Locus expansion (Zhou et al. procedure) ----

## Define loci using hit-anchored proximity expansion.
##
## A hit is a SNP with Pvalue < p_hit. All SNPs with Pvalue < p_cand within
## window_mb of any hit are retained. Loci are then split by gaps > gap_mb.
##
## @param stats      data.table with columns: chr (int), pos (int), Pvalue (num).
##                   Should already have strand-ambiguous SNPs removed.
## @param p_hit      Significance threshold defining a hit (default 1e-6).
## @param p_cand     Looser threshold for candidates (default 1e-5).
## @param window_mb  Radius around hits to include candidates (default 1.0 Mb).
## @param gap_mb     Gap size that splits loci (default 1.0 Mb).
## @param exclude_hla Drop chr6:25–34 Mb (default TRUE).
##
## @return data.table: stats rows that pass the filter, with locus_id column added.
define_loci_expanded <- function(stats, p_hit = 1e-6, p_cand = 1e-5,
                                  window_mb = 1.0, gap_mb = 1.0,
                                  exclude_hla = TRUE) {
    cands <- stats[Pvalue < p_cand]
    if (exclude_hla)
        cands <- cands[!(chr == 6L & pos >= 25e6L & pos <= 34e6L)]

    hits <- cands[Pvalue < p_hit,
                  .(chr,
                    hit_lo = pos - as.integer(window_mb * 1e6),
                    hit_hi = pos + as.integer(window_mb * 1e6))]

    if (nrow(hits) == 0L) {
        message("  define_loci_expanded: no hits found")
        return(cands[0L])
    }
    setkeyv(hits, c("chr", "hit_lo", "hit_hi"))

    cands[, pos_lo := pos]
    cands[, pos_hi := pos]
    setkeyv(cands, c("chr", "pos_lo", "pos_hi"))

    ov      <- foverlaps(cands, hits, nomatch = 0L, which = TRUE)
    kept_i  <- unique(ov$xid)
    retained <- cands[kept_i]
    cands[, c("pos_lo", "pos_hi") := NULL]

    setorder(retained, chr, pos)
    retained[, prev_pos  := shift(pos, 1L, type = "lag"), by = chr]
    retained[, new_locus := is.na(prev_pos) | (pos - prev_pos) > gap_mb * 1e6]
    retained[, locus_n   := cumsum(new_locus), by = chr]
    retained[, locus_id  := paste0("chr", chr, "_locus", locus_n)]
    retained[, c("prev_pos", "new_locus", "locus_n") := NULL]
    retained
}

## ---- Truncated eigendecomposition ----

## Returns factored form (V_k, d_inv, eff_rank) to avoid O(n^2 * eff_rank)
## materialization of the dense inverse. Apply as:
##   alpha_m = V_k %*% (d_inv * (t(V_k) %*% b))
trunc_pseudoinv <- function(S, min_eig_frac) {
    eig  <- eigen(S, symmetric = TRUE)
    tol  <- max(eig$values) * min_eig_frac
    keep <- eig$values > tol
    list(V_k      = eig$vectors[, keep, drop = FALSE],
         d_inv    = 1 / eig$values[keep],
         eff_rank = sum(keep))
}

## ---- Per-locus coefficient functions ----

## Compute LD-adjusted instrument-exposure (alpha) coefficients for one locus.
##
## @param chrom        Chromosome number (integer; non-integer returns NULL).
## @param snps         data.table with columns: rsid, Allele1, Effect, Freq1,
##                     StdErr, TotalSampleSize, pos.
##                     Allele1 must already be in the caller's reference convention.
##                     Alignment to the LD matrix A1 is done here.
## @param get_ld_fn    Closure: function(chrom, rsids) -> list(C_g, snp_df).
##                     Either wrapping get_locus_ld() or get_locus_ld_1kg().
## @param min_eig_frac Eigenvalue truncation fraction for pseudoinverse.
##
## @return Named list: chr, locus_start, locus_end, n_snps, eff_rank,
##                     alpha_hat, se.alpha_hat, sd.Z  — or NULL if no LD match.
locus_compute_alpha <- function(chrom, snps, get_ld_fn, min_eig_frac) {
    result <- get_ld_fn(chrom, snps$rsid)
    if (is.null(result[[1L]])) return(NULL)
    C_g   <- result[[1L]]
    ld_df <- as.data.table(result[[2L]])

    ld_dt  <- ld_df[, .(rsid = SNP, ld_A1 = A1, ld_A2 = A2, ld_idx = .I)]
    snps_m <- merge(ld_dt, snps, by = "rsid", sort = FALSE)
    setorder(snps_m, ld_idx)

    ## Align to LD A1: flip Effect and Freq1 where Allele1 matches LD A2
    flip <- snps_m$Allele1 == snps_m$ld_A2
    drop <- snps_m$Allele1 != snps_m$ld_A1 & snps_m$Allele1 != snps_m$ld_A2
    snps_m[flip, `:=`(Effect = -Effect, Freq1 = 1 - Freq1)]
    snps_m <- snps_m[!drop]
    if (nrow(snps_m) == 0L) return(NULL)

    keep <- snps_m$ld_idx
    C_g  <- C_g[keep, keep, drop = FALSE]
    k    <- nrow(snps_m)

    sigma_j      <- sqrt(2 * snps_m$Freq1 * (1 - snps_m$Freq1))
    alpha_u_star <- snps_m$Effect * sigma_j     # per-SD marginal effects

    if (k == 1L) {
        alpha_m  <- alpha_u_star
        eff_rank <- 1L
    } else {
        pi_out   <- trunc_pseudoinv(C_g, min_eig_frac)
        alpha_m  <- as.numeric(pi_out$V_k %*% (pi_out$d_inv * (t(pi_out$V_k) %*% alpha_u_star)))
        eff_rank <- pi_out$eff_rank
    }

    beta_m    <- alpha_m / sigma_j              # per-allele weights (eq. 7)
    alpha_hat <- sqrt(sum(beta_m^2))

    var_S <- as.numeric(t(alpha_m) %*% C_g %*% alpha_m)   # analytical Var(S)
    var_Z <- var_S / alpha_hat^2

    ## sigma_X^2 from allele freq, SE, N (eq. 9); median across SNPs
    p          <- snps_m$Freq1
    se         <- snps_m$StdErr
    eff        <- snps_m$Effect
    N          <- snps_m$TotalSampleSize
    sigma_X_sq <- median(2 * p * (1 - p) * (N * se^2 + eff^2))
    N_alpha    <- median(N)

    ## SE via Fisher information (eq. 10)
    denom    <- sigma_X_sq - var_S
    se_alpha <- if (denom <= 0) {
        warning("locus_compute_alpha: non-positive Fisher information denominator")
        NA_real_
    } else {
        1 / sqrt(N_alpha * var_Z / denom)
    }

    list(chr          = chrom,
         locus_start  = min(snps_m$pos),
         locus_end    = max(snps_m$pos),
         n_snps       = k,
         eff_rank     = eff_rank,
         alpha_hat    = alpha_hat,
         se.alpha_hat = se_alpha,
         sd.Z         = sqrt(var_Z))
}

## Compute LD-adjusted instrument-outcome (gamma) coefficients for one locus.
##
## @param chrom        Chromosome number (integer; non-integer returns NULL).
## @param snps         data.table with columns: rsid, Allele1, Effect, Freq1, pos.
##                     Allele1 in the caller's reference convention (e.g. bim a1).
## @param z_gamma      Named numeric vector of outcome GWAS z-scores (names = rsids).
##                     Must cover all rsids in snps.
## @param N_gamma      Total outcome GWAS sample size.
## @param p_case       Case fraction in outcome GWAS.
## @param get_ld_fn    Closure: function(chrom, rsids) -> list(C_g, snp_df).
## @param min_eig_frac Eigenvalue truncation fraction for pseudoinverse.
##
## @return Named list: chr, locus_start, locus_end, n_snps, eff_rank,
##                     gamma_hat, se.gamma_hat  — or NULL if no LD match.
locus_compute_gamma <- function(chrom, snps, z_gamma, N_gamma, p_case,
                                get_ld_fn, min_eig_frac) {
    result <- get_ld_fn(chrom, snps$rsid)
    if (is.null(result[[1L]])) return(NULL)
    C_g   <- result[[1L]]
    ld_df <- as.data.table(result[[2L]])

    ld_dt    <- ld_df[, .(rsid = SNP, ld_A1 = A1, ld_A2 = A2, ld_idx = .I)]
    snps_aln <- merge(ld_dt, snps, by = "rsid", sort = FALSE)
    setorder(snps_aln, ld_idx)

    ## Reorder z_gamma to LD order (z_gamma is named by rsid)
    gz <- z_gamma[snps_aln$rsid]

    ## Align to LD A1: flip both PDCD1 Effect and outcome z where Allele1 == LD A2
    flip <- snps_aln$Allele1 == snps_aln$ld_A2
    drop <- snps_aln$Allele1 != snps_aln$ld_A1 & snps_aln$Allele1 != snps_aln$ld_A2
    snps_aln[flip, `:=`(Effect = -Effect, Freq1 = 1 - Freq1)]
    gz[flip] <- -gz[flip]
    snps_aln <- snps_aln[!drop]
    gz       <- gz[!drop]
    if (length(gz) == 0L) return(NULL)

    keep <- snps_aln$ld_idx
    C_g  <- C_g[keep, keep, drop = FALSE]
    k    <- nrow(snps_aln)

    sigma_j      <- sqrt(2 * snps_aln$Freq1 * (1 - snps_aln$Freq1))
    alpha_u_star <- snps_aln$Effect * sigma_j
    gamma_u_star <- gz / sqrt(N_gamma * p_case * (1 - p_case))

    if (k == 1L) {
        alpha_m  <- alpha_u_star
        eff_rank <- 1L
    } else {
        pi_out   <- trunc_pseudoinv(C_g, min_eig_frac)
        alpha_m  <- as.numeric(pi_out$V_k %*% (pi_out$d_inv * (t(pi_out$V_k) %*% alpha_u_star)))
        eff_rank <- pi_out$eff_rank
    }

    beta_m    <- alpha_m / sigma_j
    alpha_hat <- sqrt(sum(beta_m^2))
    if (alpha_hat == 0) return(NULL)

    var_S <- as.numeric(t(alpha_m) %*% C_g %*% alpha_m)
    if (var_S <= 0) return(NULL)
    gamma_hat <- alpha_hat * sum(gamma_u_star * alpha_m) / var_S   # eq. 16
    var_Z     <- var_S / alpha_hat^2
    se_gamma  <- 1 / sqrt(N_gamma * var_Z * p_case * (1 - p_case))  # eq. 17

    list(chr          = chrom,
         locus_start  = min(snps_aln$pos),
         locus_end    = max(snps_aln$pos),
         n_snps       = k,
         eff_rank     = eff_rank,
         gamma_hat    = gamma_hat,
         se.gamma_hat = se_gamma)
}
