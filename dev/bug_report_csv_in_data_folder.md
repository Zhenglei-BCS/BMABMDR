# Bug: CSV/TXT files in `data/` are not installed with the package

## Summary

The `data/` directory contains several `.csv` and `.txt` files that are **not
accessible after package installation**. R only installs `.rda`, `.RData`, and
`.R` files from `data/`. Raw text files intended for user access (examples,
vignettes, `system.file()`) must be placed in `inst/extdata/` instead.

## Affected files

| File | Used in |
|------|---------|
| `data/EggShellThickness.csv` | Litter effect vignette (continuous clustered example) |
| `data/example_quantal.csv` | Unknown — possibly example scripts |
| `data/example_clusteredQuantal.csv` | Unknown — possibly clustered quantal examples |
| `data/examplecontinuouswithcovariate.csv` | Unknown — possibly covariate examples |
| `data/test_data.csv` | Unknown — possibly testing |
| `data/test_data_cont.txt` | Introduction vignette (individual continuous data) |

## Current behaviour

```r
system.file("data", "EggShellThickness.csv", package = "BMABMDR")
#> ""
# File does not exist in the installed package

system.file("extdata", "EggShellThickness.csv", package = "BMABMDR")
#> ""
# Not in inst/extdata/ either
```

Any code using `system.file()` or `read.csv(system.file(...))` to load these
files will silently get an empty path, then fail with "no lines available in
input" or "cannot open the connection".

## Expected behaviour

Raw data files intended for use in examples, vignettes, or by end users should
be findable via `system.file("extdata", "<filename>", package = "BMABMDR")`.

## Suggested fix

Move the files from `data/` to `inst/extdata/`:

```bash
mkdir -p inst/extdata
git mv data/EggShellThickness.csv         inst/extdata/
git mv data/example_quantal.csv           inst/extdata/
git mv data/example_clusteredQuantal.csv  inst/extdata/
git mv data/examplecontinuouswithcovariate.csv inst/extdata/
git mv data/test_data.csv                 inst/extdata/
git mv data/test_data_cont.txt            inst/extdata/
```

Then update all references to use:

```r
system.file("extdata", "EggShellThickness.csv", package = "BMABMDR")
```

Alternatively, convert each CSV to an `.rda` object and load via `data()`:

```r
EggShellThickness <- read.csv("data-raw/EggShellThickness.csv")
usethis::use_data(EggShellThickness)
# Then access via: data("EggShellThickness")
```

## R documentation reference

From [Writing R Extensions §1.1.5](https://cran.r-project.org/doc/manuals/R-exts.html#Data-in-packages):

> The `data/` directory is for data objects to be loaded via `data()`. Only
> `.rda`, `.RData`, `.tab`, `.txt`, `.csv` files are supported **but** `.csv`
> and `.txt` are loaded via `read.table()` with specific conventions. For
> arbitrary raw files, use `inst/extdata/`.

In practice, while R *can* load CSVs from `data/` via `data()`, they are subject
to strict naming and format conventions. For files accessed by `read.csv()` or
`system.file()`, `inst/extdata/` is the correct location.

## Impact

- Vignettes that reference these files fail after a clean `R CMD INSTALL`
- `system.file()` calls return `""`, leading to confusing errors
- Examples in `.Rd` files or articles that use these files will break on CRAN

## Workaround

In vignettes, use `here::here("data", "<file>")` or relative paths
(`"../data/<file>"`), which only work from the source tree — not from an
installed package.

## Environment

- BMABMDR 0.1.20
- R 4.5.1
