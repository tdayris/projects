#!/usr/bin/R

# Loading libraries
# Handle large datasets
base::library(package = dplyr, quietly = TRUE);
base::library(package = tidyr, quietly = TRUE);
base::library(package = tibble, quietly = TRUE);
base::library(package = readr, quietly = TRUE);

# Plotting behaviour
base::library(package = ggplot2, quietly = TRUE);

# Deconvolution methods
base::library(package = immunedeconv, quietly = TRUE);

# Load dataset
quanTISeq_dataset <- readr::read_tsv(
  file = snakemake@input[["fraction"]],
  colname = TRUE
);

quanTISeq_dataset <- quanTISeq_dataset %>% tidyr::gather(
  key = cell_type,
  value = value,
  2:ncol(quanTISeq_dataset)
) %>% tidyr::spread_(
  key = names(quanTISeq_dataset)[1],
  value = 'value'
);
base::print(quanTISeq_dataset %>% base::head);
base::message("Dataset loaded");

# plot histogram
grDevices::png(
  filename = snakemake@output[["histogram"]],
  width = 480,
  height = 480,
  units = "px"
);

quanTISeq_dataset %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(quanTISeq_dataset)));

dev.off();
base::message("Histogram saved");

# Dotplots
grDevices::png(
  filename = snakemake@output[["dotplots"]],
  width = 480,
  height = 480,
  units = "px"
);

quanTISeq_dataset %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x=sample, y=score, color=cell_type)) +
    geom_point(size=4) +
    facet_wrap(~cell_type, scales="free_x", ncol=3) +
    scale_color_brewer(palette="Paired", guide=FALSE) +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1));

dev.off();
base::message("Ditplots saved");
base::message("Process over");
