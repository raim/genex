---
bibliography: doc/genex.bib
---

<!-- pandoc README.md --filter pandoc-citeproc -o README.pdf -->

# Induced Gene Expression with Inducer Degradation

This git holds an implementation of the model for induced gene
expression as described in the Supporting Information `SI File 3` of
@Behle2020 (https://doi.org/10.1021/acssynbio.9b00505).

See (doc/genex.rmd) for details on the mathematical model. 
To compile in R this use:

``` R
rmarkdown::render("doc/genex.rmd", "pdf_document")
```

## Functions

The file `(R/genex.R)` holds two functions.

* `Gauss2F1` implements the described wrapper for the `hyperg_2F1`
function, avoiding convergence issues by a transformation @Forrey1997,
(https://doi.org/10.1006/jcph.1997.5794) or alternatively a solution
suggested at
https://stats.stackexchange.com/questions/33451/computation-of-hypergeometric-function-in-r
by Stephane Laurent.
* `fexpr` implements the equation for induced protein expression, accounting
for inducer half-life:

\begin{equation*}
\label{eq:vari}
y(t) = y_0 e^{-\beta t} + \frac{1- e^{-\beta t}}{\beta}\left(\ell + v\; _2\text{F}_1\left(1,\frac{\beta}{n \delta};\frac{\beta}{n \delta}+1;-e^{n\delta t}\frac{K^n}{I_0^n}\right)\right)\,.
\end{equation*}

## Examples

* `scripts/analyse.R` exemplifies effects of parameters and is sourced
for the documentation in `doc/genex.rmd`,
* `scripts/vannilate_fit.R` exemplifies how to fit a real life
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
