

# Morphological Profiling Statistics

## Introduction

This is a collection of R tools for statistical modeling morphological profiling experiments aka CellPainting.

This is built on `R` and meant to be run from within KNIME.

## Installiation

### Install R packages

    devtools::install_github("momeara/MPStats")
    install.packages("Rserve")

# Run

Start `Rserver` in R

    library(Rserve)
    Rserve(args = "--vanilla")

download data

    wget URL-TBD/ProbPos_55Plates.csv

Open `vignettes/AGM_chapter/knime_workflows/MPStats_AGM.knar.knwf` in KNIME. Set the path for the plate data "parse well scores" panel

    
# Analysis results:

## Batch effects


    