---
title: "README.md"
---


Morphological Profiling Statistics

# Introduction

This is a collection of R tools for statistical modeling morphological profiling experiments aka CellPainting.

This is built on `R` to be run from KNIME.



# Installiation

## Install R packages

    devtools::install_github("momeara/mp_stats")
    install.packages("Rserve")

# Run

Start `Rserver` in R

    library(Rserve)
    Rserve(args = "--vanilla")



    