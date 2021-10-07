# Induced Gene Expression with Inducer Degradation

This git holds an implementation of the model for induced gene
expression as described in the Supporting Information `SI File 3` of
@Behle2020 (https://doi.org/10.1021/acssynbio.9b00505).

See [doc/genex.rmd](doc/genex.rmd) for details on the mathematical
model, and how to use the model to fit protein expression data
measured in in growing cells (eg. GFP in *E. coli* in a
"platereader").  To compile in R use:

``` R
rmarkdown::render("doc/genex.rmd", "pdf_document")
```

## Usage

The file [R/genex.R](R/genex.R) holds two functions.

* `Gauss2F1` implements a wrapper for the `hyperg_2F1` function,
avoiding convergence issues by a transformation desribed by @Forrey1997
(https://doi.org/10.1006/jcph.1997.5794) or alternatively a solution
suggested at
https://stats.stackexchange.com/questions/33451/computation-of-hypergeometric-function-in-r
by Stephane Laurent.
* `fexpr` implements the equation for induced protein expression with
inducer half-life and cell growth.


## Examples

* [scripts/analyse.R](scripts/analyse.R) exemplifies effects of parameters and is sourced in [doc/genex.rmd](doc/genex.rmd),
* [scripts/vanillate_fit.R](scripts/vanillate_fit.R) exemplifies how to fit a real life
data set.


## Requirements

The R package `gsl` is required for the `hyperg_2F1` function.
This is an interface to the GNU Scientific Library. To install
on Ubuntu use:

```{bash}
## while at it you can also try to install the R package via bash:
sudo apt install libgsl-dev r-cran-gsl
```

In R:

```{R}
install.packages("gsl")
```
## 
