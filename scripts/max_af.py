#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Compute max AF in populations
"""


import pandas


# Load dataset
tsv = pandas.read_csv(
    snakemake.input['tsv'],
    sep="\t",
    header=0,
    index_col=None
)

# Add MAX AF information
tsv_iterator = zip(
    tsv["dbNSFP_ExAC_NFE_AF "],
    tsv["dbNSFP_ExAC_SAS_AF "],
    tsv["dbNSFP_ExAC_Adj_AF "],
    tsv["dbNSFP_ExAC_AFR_AF "],
    tsv["dbNSFP_ExAC_AF "],
    tsv["dbNSFP_ExAC_FIN_AF "],
    tsv["dbNSFP_1000Gp3_EUR_AF "],
    tsv["dbNSFP_ExAC_AMR_AF "],
    tsv["dbNSFP_ExAC_EAS_AF "],
    tsv["dbNSFP_ESP6500_EA_AF "],
    tsv["dbNSFP_1000Gp3_AMR_AF "],
    tsv["dbNSFP_1000Gp3_AF "],
    tsv["dbNSFP_1000Gp3_EAS_AF "],
    tsv["dbNSFP_1000Gp3_EUR_AF "],
    tsv["dbNSFP_1000Gp3_AFR_AF "],
    tsv["dbNSFP_1000Gp3_SAS_AF"]
)

tsv["MAX_AF"] = [max(i) for i in tsv_iterator]

# Save results
tsv.to_csv(snakemake.output["tsv"], sep="\t", index=False)
