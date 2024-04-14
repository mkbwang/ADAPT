# ADAPT

## Installation

To install the package with the vignette generated, type in R

```{r, eval=FALSE}
install.packages("devtools")
library(devtools)
install_github("mkbwang/ADAPT", build_vignettes=TRUE)
```

_The authors are preparing to submit this package to bioconductor._


## Introduction

`ADAPT` stands for "Analysis of microbiome differential abundance by pooling Tobit models". There are two main innovations:

* `ADAPT` treats zero counts as left censored observations and use Tobit models to model the log count ratios.
* `ADAPT` has an innovative way of selectin non-differentially abundant taxa as reference taxa. It then uses the reference taxa to find the differentially abundant ones.

To learn more about the usage of `ADAPT`, type in R console

```{r, eval=FALSE}
browseVignettes(package="ADAPT")
```

