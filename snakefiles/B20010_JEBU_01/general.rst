Material and methods
====================

Material and Methods
====================

The raw reads' quality were assessed with FastQC, and let us conclude no further trimming was needed. The quantification process was performed on Ensembl GRCh38.86 cdna, with Salmon. While a perfect hash was used while indexing, and duplicate targets were kept. The quantification process used 100 bootstrap rounds, mapping validation, and both sequence bias and GC bias taken into account. Transcript to gene count aggregation was performed by Salmon using the corresponding GTF from ensembl.

Count tables were aggregated and control plots were made with in-house scripts.

The whole pipeline was powered by Snakemake and the Snakemake wrappers to ensure a fair distribution, repeatability, reproduction and understanding. These scripts are available on demand.


Citations
---------

This pipeline is available on github.com/tdayris/projects/Snakefile/B20010_JEBU_01/

FastQC
  ANDREWS, Simon, et al. FastQC: a quality control tool for high throughput sequence data. 2010.

  Why FastQC? FastQC is a quite popular tool among the field of bioinformatics and genomics. Cited more than 4.000 times since its publication, this tool performs reliable quality controls on sequenced reads.

  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Salmon
  Patro, Rob, et al. “Salmon provides fast and bias-aware quantification of transcript expression.” Nature Methods (2017). Advanced Online Publication. doi: 10.1038/nmeth.4197.

  Why Salmon? Cited more than 1.200 the last three years, Salmon is a very powerful tool to quantify transcripts taking into account a wide range of bias (either biological, or bioinformatical ones).

  https://salmon.readthedocs.io/

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
