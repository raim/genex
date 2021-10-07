---
bibliography: doc/genex.bib
---

<!-- pandoc README.md --filter pandoc-citeproc -o README.pdf -->

# Induced Gene Expression with Inducer Degradation

This implements the model for induced gene expression as
described in Behle et al. 2020 (https://doi.org/10.1021/acssynbio.9b00505).

See (doc/genex.rmd) for details on the mathematical model (compile in Rstudio
or with `Rscript -e 'library(rmarkdown); rmarkdown::render("genex.rmd", "pdf_document")'`).

The file `(R/genex.R)` implements the described wrapper for the `hyperg_2F1` function,
avoiding convergence issues by a transformation @Forrey1997,
and implements the model outlined in 

# Usage

## Requirements

The R package `gsl` is required for the `hyperg_2F1` function.
This is an interface to the GNU Scientific Library. To install
on Ubuntu use:

```{bash}
sudo apt install libgsl-dev
## while at it you can try to install the R package via bash:
#sudo apt install -y r-cran-gsl
```

In R:

```{R}
install.packages("gsl")
```
## 
