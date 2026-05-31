# Scaling bug in compute_alpha.R and compute_gamma.R

## Diagnosis (chr6_locus2, single SNP rs72928038)

| Quantity | Value |
|----------|-------|
| Effect (per A allele, from GWAS) | 0.041732 SD_X |
| StdErr | 0.007960 |
| N | 48918 |
| Freq1 (A allele) | 0.1786 |
| SD_g from 1000G ref panel | 0.538528 |
| sqrt(2*f*(1-f)) from Freq1 | 0.541668 |

### Step-by-step

1. **Raw coefficient (alpha_u):** Effect = 0.041732 SD_X per A allele  
2. **LD-adjusted (alpha_m):** k=1, so alpha_m = alpha_u = 0.041732 (unchanged)  
3. **Reference panel scores:**  
   - SD_g = 0.538528 (SD of raw genotype in 1000G EUR)  
   - Z = G_std (unit variance in ref panel, Var(Z) = 1.000)  
4. **alpha_hat computed by code:** |Effect| = **0.041732**  
   **alpha_hat correct:** |Effect| × SD_g = **0.022474**  
   - Code uses per-allele effect; Z is unit-variance → coefficient should be per-unit-Z = per-SD_g

### SE comparison

| Method | se_alpha |
|--------|----------|
| Code (uses alpha_hat = \|Effect\|) | 0.004309 |
| Corrected (alpha_hat = \|Effect\| × SD_g) | 0.004312 |
| StdErr × SD_g | 0.004287 |
| StdErr (raw GWAS) | 0.007960 |

The SE formula gives ≈ StdErr × SD_g regardless (because Effect² ≪ sigma_X²), but alpha_hat is wrong by a factor of 1/SD_g ≈ 1.86.

### Consequence for z-score significance of alpha_hat

| Version | alpha_hat | se_alpha | z |
|---------|-----------|----------|---|
| Code | 0.041732 | 0.004309 | **9.68** (wrong) |
| Correct | 0.022474 | 0.004312 | **5.21** (correct) |
| GWAS z = Effect/StdErr | | | 5.24 ✓ |

The code inflates the z-score for alpha_hat by ~√(1/(2f(1-f))) ≈ 1.85.  
The ratio theta_IV = gamma_hat / alpha_hat is **unaffected** because gamma_hat has the same scaling error.

---

## Root cause

The code uses `alpha_u = snps$Effect` (per-allele effect) but the pseudoinverse uses
`Sigma_g = cor(G_std)` (the correlation matrix of *standardized* genotypes). The correct
`alpha_u` to pair with Sigma_g (correlation matrix) is the *per-SD-genotype* effect:

```r
alpha_u_correct <- snps$Effect * col_sds   # col_sds = SD of G in ref panel
```

Similarly in compute_gamma.R, `gamma_u` is computed as log-OR per raw allele:
```r
# current (per-allele):
se_logOR <- 1 / sqrt(N_gamma * 2 * f * (1 - f) * p_case * (1 - p_case))
gamma_u  <- z_t1d * se_logOR

# correct (per-SD genotype):
se_logOR <- 1 / sqrt(N_gamma * p_case * (1 - p_case))
gamma_u  <- z_t1d * se_logOR
```

---

## Fix required in compute_alpha.R

In `process_locus()`, replace the k=1 and k>1 branches:

```r
## Current (wrong):
if (k == 1L) {
    alpha_m  <- snps$Effect
    eff_rank <- 1L
} else {
    ...
    alpha_m  <- as.numeric(pi_out$inv %*% snps$Effect)
}

## Corrected:
if (k == 1L) {
    alpha_m  <- snps$Effect * col_sds
    eff_rank <- 1L
} else {
    ...
    alpha_m  <- as.numeric(pi_out$inv %*% (snps$Effect * col_sds))
}
```

Expected result for chr6_locus2 after fix:
- alpha_hat = 0.022474 (vs 0.041732 currently)
- se_alpha  ≈ 0.004312 (essentially unchanged)
- z ≈ 5.21 (correct; matches GWAS z ≈ 5.24)

## Fix required in compute_gamma.R

In `process_locus()`, change the gamma_u conversion:

```r
## Current (per-allele):
se_logOR <- 1 / sqrt(N_gamma * 2 * f * (1 - f) * p_case * (1 - p_case))
gamma_u  <- gamma_snps$z_t1d * se_logOR

## Corrected (per-SD genotype):
se_logOR <- 1 / sqrt(N_gamma * p_case * (1 - p_case))
gamma_u  <- gamma_snps$z_t1d * se_logOR
```

Note: the alpha_u fix in compute_gamma.R needs the same col_sds multiplication:
```r
alpha_u  <- alpha_snps$Effect * col_sds   # was: alpha_snps$Effect
```

---

## Impact on MR-Hevo results

- theta = gamma_hat / alpha_hat: **unchanged** (both off by same SD_g factor)
- Individual alpha_hat values: **wrong by factor 1/SD_g** (~1.7–2x for common variants)  
- Significance of alpha_hat: **inflated** (~2x z-score)  
- Reported SE: **correct** (gives StdErr × SD_g in both versions)  
- MR-Hevo posterior for theta: **not affected** (uses ratio)
