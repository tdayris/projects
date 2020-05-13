### Material and methods

Fastq_ files quiality were assessed with FastQC_, showing good overall quality. Raw reads were mapped against GRCm38_, using BWA_. They were sorted, corrected, realigned and indexed using Picard_, GATK_ and Samtools according to the pipeline wes-mapping-bwa-gatk_.

Calling was performed with BCFTools, Mustect2, and HaplotypeCaller. Consensus calls were gathered and filtered with BCFTools_. The annotation was performed with SnpEff_, and the format conversion from VCF_ to TSV_ was done with SnpSift_.

### Results

The final TSV_ tables contain the following columns:

+-------------------+--------------------------------------------+
| Column            | Content                                    |
+===================+============================================+
| CHROM             | Chromosome (or un-assembled contig)        |
+-------------------+--------------------------------------------+
| POS               | Position within the Chromosome             |
+-------------------+--------------------------------------------+
| REF               | Reference allele on the sequence           |
+-------------------+--------------------------------------------+
| ALT               | Alternative(s) alleles                     |
+-------------------+--------------------------------------------+
| GT                | Genotype                                   |
+-------------------+--------------------------------------------+
| ANN[*].GENE       | Gene name                                  |
+-------------------+--------------------------------------------+
| ANN[*].EFFECT     | Effect of the mutation                     |
+-------------------+--------------------------------------------+
| ANN[*].IMPACT     | Mutation impact prediction                 |
+-------------------+--------------------------------------------+
| ANN[*].BIOTYPE    | Transcript biotype                         |
+-------------------+--------------------------------------------+
| ANN[*].HGVS_C     | c. notation                                |
+-------------------+--------------------------------------------+
| ANN[*].HGVS_P     | p. notation                                |
+-------------------+--------------------------------------------+
| ANN[*].CDNA_POS   | Position in the cDNA                       |
+-------------------+--------------------------------------------+
| DP                | Read depth at the given position           |
+-------------------+--------------------------------------------+
| ANN[*].GENEID     | Gene ID in the databases                   |
+-------------------+--------------------------------------------+


.. _Fastq: https://en.wikipedia.org/wiki/FASTQ_format
.. _FastQC: Andrews, Simon. "FastQC: a quality control tool for high throughput sequence data." (2010). Cited more than 4000 times.
.. _GRCm38: Zerbino, Daniel R., et al. "Ensembl 2018." Nucleic acids research 46.D1 (2018): D754-D761. Cited more than 1000 time per year.
.. _BWA: http://bio-bwa.sourceforge.net/bwa.shtml#13. Cited more than 10000 times.
.. _Picard: DePristo, Mark A., et al. "A framework for variation discovery and genotyping using next-generation DNA sequencing data." Nature genetics 43.5 (2011): 491. Cited more than 7000 times.
.. _GATK: McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303. Cited more than 12000 times.
.. _Samtools: Li, Heng, et al. "The sequence alignment/map format and SAMtools." Bioinformatics 25.16 (2009): 2078-2079. Cited more than 24000 times.
.. _wes-mapping-bwa-gatk: https://github.com/tdayris-perso/wes-mapping-bwa-gatk
.. _BCFTools:
.. _Mutect2: do Valle, Ítalo Faria, et al. "Optimized pipeline of MuTect and GATK tools to improve the detection of somatic single nucleotide polymorphisms in whole-exome sequencing data." BMC bioinformatics 17.12 (2016): 341. Cited more than 1000 times.
.. _HaplotypeCaller: Poplin, Ryan, et al. "Scaling accurate genetic variant discovery to tens of thousands of samples." BioRxiv (2017): 201178. Cited more than 1000 times.
.. _SnpEff: Cingolani, Pablo, et al. "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3." Fly 6.2 (2012): 80-92. Cited more than 4000 times.
.. _VCF: https://en.wikipedia.org/wiki/Variant_Call_Format
.. _TSV: https://en.wikipedia.org/wiki/Tab-separated_values
.. _SnpSift: Ruden, Douglas Mark, et al. "Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift." Frontiers in genetics 3 (2012): 35. Cited more than 400 times.
