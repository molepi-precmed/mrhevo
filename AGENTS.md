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
