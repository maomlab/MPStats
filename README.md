![MPStats Logo](MPStats_logo.png "Logo")

# Morphological Profiling Statistics

## Introduction

This is a collection of R tools for statistical modeling morphological profiling experiments aka CellPainting.

This is built on `R` and meant to be run from within KNIME.

## Installiation

### Install R packages

    devtools::install_github("momeara/MPStats")
    install.packages("Rserve")

# Run AGM Chapter Vignette as a batch

    git clone git@github.com:momeara/MPStats.git
    cd MPStats/vignettes/AGM_chapter
    make clean
    make run_analysis
    
# Run AGM Chapter Vignette from within KNIME

Start `Rserver` in R

    library(Rserve)
    Rserve(args = "--vanilla")

1.  Open `vignettes/AGM_chapter/knime_workflows/MPStats_AGM.knar.knwf` in KNIME.
2.  Set the path for the plate data "parse well scores" panel
3.  Click `Excute All` from the commandbar


    
