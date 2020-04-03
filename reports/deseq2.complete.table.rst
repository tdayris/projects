This is the results of DESeq2, it contains the differential analysis results. With this file, you will be able to see either differentially expressed genes, and the others: this table contains both of them!

This differential analysis has been done with the following statistical formulae: `{{snakemake.wildcards.formulae}}`

It contains the following columns:

+----------------+-----------------------------------------------+
| Column         | Content                                       |
+================+===============================================+
|                | The gene identifier within Ensembl database   |
+----------------+-----------------------------------------------+
| baseMean       | The mean count for this gene                  |
+----------------+-----------------------------------------------+
| log2FoldChange | The log2(MeanExpression_C1/MeanExpression_C2) |
+----------------+-----------------------------------------------+
| lfcSE          | The shrinked value of the log2(FC)            |
+----------------+-----------------------------------------------+
| stat           | The stat value used for the linear model      |
+----------------+-----------------------------------------------+
| pvalue         | The raw P-value (should not be considered)    |
+----------------+-----------------------------------------------+
| padj           | The Adjusted P-Value                          |
+----------------+-----------------------------------------------+

The values are alphabetically sorted according to the name of the gene.

This file is a TSV file (Tab Separated Values). It can be opened by all tabular editor such as Excel, LibreOffice Calc, etc. Thy can do it, even if your computer tell you otherwise. Open your favorite tabular editor, then click and drag the TSV file into the opened tabular editor.

This file is quite big. You may try to open it, but it shall slow down your computer. If you want to see expression values between a set of genes, please ask you bioinformatician : we can draw box plots of any genes you're interested in, for any number of genes you might be interested in (1, 100, 1000, all).
