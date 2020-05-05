#!/usr/bin/R

# This script performs quantiseq deconvolution

base::library(package="readr", quietly=TRUE);
base::library(package="immunedeconv", quietly=TRUE);


# Load counts
quant <- readr::read_tsv(
    file=snakemake@input[["quant"]],
    col_names=TRUE
);
base::message("Libraries and counts loaded");


# Run deconvolution
quantiseq <- immunedeconv::deconvolute(
  quant,
  "quantiseq",
  tumor=snakemake@params[["tumor"]]
);
base::message("Deconvolution performed");


# Save results
base::saveRDS(
  quantiseq,
  file=snakemake@output[["rds"]]
);
base::message("Process over");
