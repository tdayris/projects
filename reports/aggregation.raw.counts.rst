This file contains rw counts. raw counts are *not* normalized counts : you *cannot* compare counts. Neither between multiple samples, nor between genes within a single sample.

Thus, you cannot say a gene is diferentially expressed based on this file: this is done by another analysis. Any variation visible among the count values may be misleading, only an adjusted P-Value can help you conclude (considering a given risk) if a gene is differentially expressed or not.

This file is a TSV file (Tab Separated Values). It can be opened by all tabular editor such as Excel, LibreOffice Calc, etc. Thy can do it, even if your computer tell you otherwise. Open your favorite tabular editor, then click and drag the TSV file into the opened tabular editor.

This file is quite big. You may try to open it, but it shall slow down your computer. If you want to see expression values between a set of genes, please ask you bioinformatician : we can draw box plots of any genes you're interested in, for any number of genes you might be interested in (1, 100, 1000, all).
