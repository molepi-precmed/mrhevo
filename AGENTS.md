# AGENTS.md - Guidelines for MR-Hevo Package Development

This file provides guidelines for agentic coding agents working on the MR-Hevo R package.

## 1. Build/Lint/Test Commands

### Package Check
```bash
# Full R CMD check
R CMD check --no-manual --no-vignettes /path/to/mrhevo

# Quick syntax check
R --vanilla -e "parse(file='R/mrhevo.R')"
R --vanilla -e "parse(file='R/functions.mrhevo.R')"
```

### Running Tests
```bash
# Run all tests
R CMD check --no-manual

# Run a single test file
R --vanilla -e "devtools::test_file('tests/testthat/test-functions.R')"

# Run specific test
R --vanilla -e "testthat::test_file('tests/testthat/test-functions.R', 'test_name')"
```

### Installation
```bash
# Install dependencies
R --vanilla -e "install.packages(c('data.table', 'rstan', 'ggplot2', 'car', 'crayon', 'cowplot', 'ggrepel', 'bayesplot'))"

# Load package locally
R --vanilla -e "devtools::load_all('.')"
```

## 2. Code Style Guidelines

### File Organization
- Package documentation in `R/mrhevo.R` (package-level doc, hooks)
- All exported functions in `R/functions.mrhevo.R`
- Keep Roxygen2 documentation inline with functions
- Stan models in root directory: `MRHevo_*.stan`

### Imports
- Use roxygen2 `@import` tags for package imports in NAMESPACE
- Use `@importFrom pkg func` for specific functions
- Add all dependencies to DESCRIPTION Imports field
- Available packages: crayon, data.table, ggplot2, ggrepel, car, cowplot, bayesplot, rstan

### Naming Conventions
- **Functions**: snake_case (e.g., `get_coeffratios`, `run_mrhevo.sstats`)
- **Arguments**: snake_case (e.g., `alpha_hat`, `se.alpha_hat`)
- **Variables**: snake_case (e.g., `theta_IV`, `inv.var`)
- **Constants**: UPPERCASE (e.g., `J` for number of instruments)
- **Private/helper functions**: prefix with dot (e.g., `.onAttach`)

### Code Formatting
```r
# Good
function_name <- function(arg1, arg2,
                         arg3 = default) {
    result <- compute(arg1)
    return(result)
}

# Bad - inconsistent spacing
function_name<-function(arg1,arg2){
result<-compute(arg1)
return(result)}
```

### Indentation
- Use 4 spaces for indentation
- Align function arguments across lines
- Break long lines at ~80 characters

### Roxygen2 Documentation
All exported functions must have:
- Description (first line)
- @param for each argument
- @return for return value
- @export tag
- @import tags for package dependencies

```r
#' Function title with brief description.
#'
#' Longer description if needed.
#'
#' @param arg1 Description of arg1.
#' @param arg2 Description of arg2.
#'
#' @return Description of return value.
#'
#' @import ggplot2 data.table
#' @export
function_name <- function(arg1, arg2) {
    # code
}
```

### Error Handling
- Use `stopifnot()` for argument validation
- Provide informative error messages
- Use `warning()` for non-fatal issues

```r
# Good
stopifnot(length(alpha_hat) == length(gamma_hat))
stopifnot(fraction_pleio >= 0.01 & fraction_pleio <= 0.99)

# With custom message
if (length(alpha_hat) != length(gamma_hat)) {
    stop("alpha_hat and gamma_hat must have same length")
}
```

### Data Structures
- Use `data.table` for tabular data
- Prefer `data.table()` over `data.frame()`
- Use `setDT()` for conversion
- Reference columns without quotes in data.table syntax

### Statistical Functions
- Use `matrixStats` for efficient matrix operations
- Use `rstan` for Bayesian inference (Stan models)
- Return results as `data.table` for tabular output

### Plotting
- Use `ggplot2` for all plots
- Create separate exportable plotting functions
- Always use `theme_bw()` or `theme_minimal()`
- Label axes clearly with units
- Use `ggrepel` for point labels when needed

### Testing (testthat)
- Place tests in `tests/testthat/`
- Name test files as `test-*.R`
- Use `expect_*` assertions
- Test both success and error cases

```r
test_that("function returns correct structure", {
    result <- my_function(input)
    expect_type(result, "list")
    expect_named(result, c("estimate", "se"))
})
```

## 3. Common Patterns

### Stan Model Execution
```r
rstan_options(auto_write = TRUE)
stan_model(file = file.path(model.dir, "model.stan"), verbose = FALSE)
rstan::sampling(object = stanmodel, data = data.stan, ...)
```

### Parallel Processing
```r
# Register doParallel (done in .onAttach)
options(cores = min(floor(parallel::detectCores() / 2), 10))
doParallel::registerDoParallel()

# Use foreach for parallel loops
foreach(i = 1:n, .combine = rbind) %dopar% {
    # computation
}
```

### Package Development Workflow
1. Make changes to R code in `R/`
2. Test syntax: `R --vanilla -e "parse(file='R/functions.mrhevo.R')"`
3. Run tests: `R CMD check`
4. Commit and push to GitHub

## 4. Important Notes

- The package uses RStan for Bayesian inference ( Stan models required)
- Stan model files must be in the package root or specified via `model.dir`
- The package targets R >= 4.0
- Use roxygen2 to regenerate NAMESPACE after adding exports

## 5. Computing instrument-exposure coefficients from summary statistics

The method is described in `theorymethods.Rmd` (rendered as `theorymethods.pdf`), section
"Constructing scalar instruments from multiple SNPs".  The script `compute_alpha.R`
implements task (1): given GWAS summary statistics for SNP-exposure associations, compute
for each locus the scalar instrument-exposure coefficient `alpha_hat` and its standard error.

### Inputs

| File | Description |
|------|-------------|
| `pdcd1_stats/OID00791_OID21396_SCALLOP_UKBB_MA_rsidannotated_filtered.tsv` | Meta-analysis summary statistics for PDCD1 protein (SCALLOP + UKBB); ~13k SNPs, ~47–50k N per variant, hg38 |
| `refpop/kg.2020.hg38.eur.bed/.bim/.fam` | 1000 Genomes EUR reference panel; 503 samples, 33M SNPs, hg38 |

Two meta-analyses are available (OID00791 and OID01098); OID00791 is used because it has
consistent sample sizes throughout.  OID01098 has some variants with N as low as 860.

### Algorithm

1. **Load summary statistics** (`chr`, `pos`, `rsid`, `Allele1`, `Allele2`, `Freq1`,
   `Effect`, `StdErr`, `TotalSampleSize`).

2. **Define loci**: sort by chr/pos; a gap > `gap_mb` (default 1 Mb) between adjacent
   SNPs on the same chromosome starts a new locus.

3. **Match to reference panel**: read the `.bim` file; join on `rsid`; align alleles
   (flip `Effect` and `Freq1` when the effect allele matches `a2` in the bim).
   Drop strand-ambiguous SNPs (A/T, C/G) and SNPs absent from the reference panel.

4. **Open reference panel** via `BEDMatrix` (memory-mapped; does not load the full
   4 GB `.bed` file into RAM).

5. **Per-locus computation**:
   - Extract genotype matrix *G* [503 × k]; impute missing values to column means;
     standardise columns (zero mean, unit SD).
   - Compute correlation matrix **Σ_g** = `cor(G_std)`.
   - Compute truncated pseudoinverse of **Σ_g** via eigendecomposition, retaining only
     eigenvalues above `min_eig_frac` (default 0.01) × max eigenvalue.  This handles
     ill-conditioned matrices such as the HLA region on chromosome 6, where many SNPs
     are in high LD.  The number of retained eigenvalues (`eff_rank`) is reported for
     each locus.
   - Compute multivariable coefficients **α_m** = **Σ_g**⁻¹ **α_u**.
   - `alpha_hat` = ‖**α_m**‖.
   - Estimate Var(*Z*) empirically from reference panel scores *Z_i* = (*G_i* · **α_m**) / ‖**α_m**‖.
   - Estimate *σ_X*² per SNP from allele frequency, SE, and N (formula in `theorymethods.Rmd`);
     take the median across the locus.
   - Compute SE via Fisher information:
     `se.alpha_hat = 1 / sqrt(N * Var(Z) / (sigma_X^2 - alpha_hat^2 * Var(Z)))`.

6. **Output**: `alpha_pdcd1_OID00791.rds` — a `data.table` with one row per locus and
   columns `qtlname`, `chr`, `locus_start`, `locus_end`, `n_snps`, `eff_rank`,
   `alpha_hat`, `se.alpha_hat`.  This feeds into `coeffs.dt` for the MR-Hevo
   Bayesian analysis (see `runmrhevo.R`).

### Key parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `gap_mb` | 1.0 | Mb gap between SNPs that defines a new locus boundary |
| `min_eig_frac` | 0.01 | Eigenvalue truncation threshold for pseudoinverse |

### Running

```bash
cd /mnt/pmd/mrhevo
Rscript compute_alpha.R
```

Requires R packages `data.table` and `BEDMatrix` (`install.packages("BEDMatrix")`).
Reading the 33M-row bim index takes ~30 s; per-locus processing is then fast except
for large, high-LD loci (e.g. HLA).
