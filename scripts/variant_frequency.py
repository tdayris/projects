#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Compute variant frquency based on a TSV file with given POS, CHR, and ALT
"""

import pandas
import collections

# Load dataset
tsv = pandas.read_csv(
    snakemake.input['tsv'],
    sep="\t",
    header=0,
    index_col=None
)


# Create variants identifier
tsv_iterator = zip(
    tsv[snakemake.params["chr"]],
    tsv[snakemake.params["pos"]],
    tsv[snakemake.params["alt"]]
)

tsv["variant_key"] = [
    f"{chr}_{pos}_{alt}"
    for chr, pos, alt in tsv_iterator
]

# Count variant occurences
counts = collections.Counter(
    tsv["variant_key"]
)

tsv["Recurence"] = [
    counts["key"]
    for key in tsv["variant_key"]
]

# Remove identifier
del tsv["variant_key"]

# Reorder columns
cols = tsv.columns.tolist()[:-1]
cols = cols[0] + ["Recurence"] + cols[1:]
tsv = tsv[cols]


# Save results
tsv.to_csv(
    snakemake.output["tsv"],
    sep="\t",
    index=False
)
