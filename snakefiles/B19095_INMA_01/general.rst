Material and Methods
====================

The raw reads' quality were assessed with FastQC, and let us conclude no further trimming was needed. The quantification process was performed on Ensembl GRCh38.86 cdna, with Salmon. While a perfect hash was used while indexing, and duplicate targets were kept. The quantification process used 100 bootstrap rounds, mapping validation, and both sequence bias and GC bias taken into account. Count aggregation was performed by tximport, with default salmon parameters, inferential replicates importation and transcript to gene counts aggregation. After a round of PCA analysis with pcaExplorer, the differential analysis was performed with DESeq2, using local fit regression.

The raw reads were also mapped on Ensembl GRCh38 with STAR. Encode specific parameters were used, and basic two-pass mapping. Quality controls were performed on these BAM files by Samtools.

T and B cell repertoire analysis from RNA-Seq was performed with MiXCR, following the non-targeted data recommendations: two rounds of partial assembly and final extension.

Quality controls were gathered by MultiQC for an efficient and quicker overview.

The whole pipeline was powered by Snakemake and the Snakemake wrappers to ensure a fair distribution, repeatability, reproduction and understanding. These scripts are available on demand.

Citations
---------

This pipeline is available on github.com/tdayris/projects/Snakefile/B19095_INMA_01/

FastQC
  ANDREWS, Simon, et al. FastQC: a quality control tool for high throughput sequence data. 2010.

  Why FastQC? FastQC is a quite popular tool among the field of bioinformatics and genomics. Cited more than 4.000 times since its publication, this tool performs reliable quality controls on sequenced reads.

  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Salmon
  Patro, Rob, et al. “Salmon provides fast and bias-aware quantification of transcript expression.” Nature Methods (2017). Advanced Online Publication. doi: 10.1038/nmeth.4197.

  Why Salmon? Cited more than 1.200 the last three years, Salmon is a very powerful tool to quantify transcripts taking into account a wide range of bias (either biological, or bioinformatical ones).

  https://salmon.readthedocs.io/

tximport
  Love, Michael I., Charlotte Soneson, and Mark D. Robinson. "Importing transcript abundance datasets with tximport." dim (txi. inf. rep $ infReps $ sample1) 1.178136 (2017): 5.

  Why tximport? While this tool never received proper publication, it has been made by the team of DESeq2 to ease the use of their tool with Salmon itself.

  https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

DESeq2
  Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome biology 15.12 (2014): 550.

  Why DESeq2? Cited more than 15.000 times, this tool is very well known and very often tested among both other bioinformatics tools, and biology qRT-PCRs.

  https://bioconductor.org/packages/release/bioc/html/DESeq2.html


pcaExplorer
  Marini, Federico, and Harald Binder. "pcaExplorer: an R/Bioconductor package for interacting with RNA-seq principal components." BMC bioinformatics 20.1 (2019): 331.

  Quite new, this tool is here to help assessing right factors to possible statistical forumlas. By itself, this tool does not perform any computation.

  https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html

MiXCR
  Bolotin, Dmitriy A., et al. "MiXCR: software for comprehensive adaptive immunity profiling." Nature methods 12.5 (2015): 380.

  Cited about 500 times, this tool performs immunity profiling on RNA-Seq data. It handles fastq file from a dedicated mapping step up to the clones assembly.

  https://mixcr.readthedocs.io/en/master/

STAR
  Dobin, Alexander, et al. "STAR: ultrafast universal RNA-seq aligner." Bioinformatics 29.1 (2013): 15-21.

  Why STAR? Cited more than 10.500 times, STAR is a very efficient gap-aware read mapper. While using a large amount of informatical resources, it performs very well in terms of time and results quality.

  https://github.com/alexdobin/STAR

Samtools
  Li, Heng, et al. "The sequence alignment/map format and SAMtools." Bioinformatics 25.16 (2009): 2078-2079.

  Why Samtools? Cited more than 24.000 times, samtools is a staple of the bioinformatician cuisine.

  http://samtools.sourceforge.net/

MultiQC
  EWELS, Philip, MAGNUSSON, Måns, LUNDIN, Sverker, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016, vol. 32, no 19, p. 3047-3048.

  Why MultiQC? MultiQC is a very efficient tool when it comes to quality gathering. It has been cited more than 500 times in a very wide range of journals including Nature, Bioinformatics, Cell, etc.

  https://multiqc.info/

Snakemake
  Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

  Why Snakemake? Snakemake is a very popular workflow manager in data science and bioinformatics. It has about three new citations per week within the scopes of biology, medicine and bioinformatics.

  https://snakemake.readthedocs.io/
  https://snakemake-wrappers.readthedocs.io/

More information
----------------

The whole pipeline is displayed below. At the bottom of this page, you can find run times, creation dates, and all the command lines, scripts and environments used to analyze this project.
