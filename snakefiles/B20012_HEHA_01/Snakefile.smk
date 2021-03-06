import pandas

from pathlib import Path
from typing import Any, Dict, List

design_path = Path("design.csv")
swv="0.50.4"
git = "https://raw.githubusercontent.com/tdayris-perso/snakemake-wrappers"

design = pandas.read_csv(
    design_path,
    sep="\t",
    header=0,
    index_col=0
)

design = design[design["Upstream_file"] != "None"]

design["R1_Copy"] = [
    f"raw_data/{Path(f).name}" for f in design["Upstream_file"]
]

design["R2_Copy"] = [
    f"raw_data/{Path(f).name}" for f in design["Downstream_file"]
]

design_dict = design.to_dict()

report: "general.rst"

def sample_pair_w(wildcards: Any) -> Dict[str, str]:
    """
    Given a sample name, this function returns a pair of fastq files
    """
    return {
        "fq1": design_dict["R1_Copy"][wildcards.sample],
        "fq2": design_dict["R2_Copy"][wildcards.sample]
    }


def salmon_sample_pair_w(wildcards: Any) -> Dict[str, str]:
    """
    Same as above, rename keys
    """
    return {
        "r1": design_dict["R1_Copy"][wildcards.sample],
        "r2": design_dict["R2_Copy"][wildcards.sample]
    }


samplsdict = {
    os.path.basename(f): f
    for f in chain(
        design["R1_Copy"],
        design["R2_Copy"]
    )
}

copy_fq_dict = {
    os.path.basename(f): f
    for f in chain(
        design["Upstream_file"],
        design["Downstream_file"]
    )
}

localrules: download_fasta, download_fasta_cdna, download_gtf, copy_fq

rule all:
    input:
        multiqc = "qc/multiqc.html",
        tximport = "deseq/tximport.RDS",
        pairwise_scatterplot = "figures/pairwise_scatterplot.png",
        count_table = "salmon/aggregated/TPM.counts.tsv",
        boxplot = "figures/box_counts.png",
        pca = expand(
            "figures/PCA/PCA_{factor}_PC1_PC2.png",
            factor = [
                "Localization", "Gender", "Early_metastasis", "Treatment",
                "KRAS_Mutation", "RAF_Mutation", "NRAS_Mutation",
                "TP53_Mutation", "BRAF_Mutation", "APC_Mutation",
                "PI3KCA_Mutation", "RAS_Mutation", "FAM_Mutation"
            ]
        ),
        deseq2_table = expand(
            "deseq2/TSV/DESeq2_{factor}",
            factor=["Liver", "Peritoneum"]
        )
    message:
        "Finishing pipeline"

rule copy_fq:
    input:
        lambda wildcards: copy_fq_dict[wildcards.fq]
    output:
        "raw_data/{fq}"
    message:
        "Copying {wildcards.fq}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(128 * attempt, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(30 * attempt, 120)
        )
    wildcard_constraints:
        fq = r"[^/]+"
    conda:
        "../envs/bash.yaml"
    log:
        "logs/copy/{fq}.log"
    shell:
        "cp -v {input} {output} > {log} 2>&1"


rule fastqc:
    input:
        lambda wildcards: samplsdict[wildcards.sample]
    output:
        html = "qc/fastqc/{sample}.html",
        zip = "qc/fastqc/{sample}_fastqc.zip"
    message:
        "Controlling quality of {wildcards.sample}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(1024 + 1024 * attempt, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(30 * attempt, 120)
        )
    params:
        ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        f"{swv}/bio/fastqc"


rule download_gtf:
    output:
        temp("resources/Homo_sapiens.gtf")
    message:
        "Downloading GTF from ensembl"
    params:
        species = "homo_sapiens",
        release = "98",
        build = "GRCh38",
        fmt = "gtf"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(128 * attempt, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(25 * attempt, 120)
        )
    log:
        "logs/ensembl_annotation/get_gtf.log"
    wrapper:
        f"{swv}/bio/reference/ensembl-annotation"


rule download_fasta:
    output:
        temp("resources/Homo_sapiens.fasta")
    message:
        "Downloading Fasta from ensembl"
    params:
        species = "homo_sapiens",
        release = "98",
        datatype = "dna",
        build = "GRCh38"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(128 * attempt, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(25 * attempt, 120)
        )
    log:
        "logs/ensembl_annotation/get_genome.log"
    wrapper:
        f"{swv}/bio/reference/ensembl-sequence"


rule star_index:
    input:
        fasta = "resources/Homo_sapiens.fasta",
        gtf = "resources/Homo_sapiens.gtf"
    output:
        directory("star/index/Homo_sapiens")
    message:
        "Indexing GRCh38 with STAR"
    threads:
        10
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(35840 + 10240 * attempt, 51200)
        ),
        time_min = (
            lambda wildcards, attempt: min(35 * attempt, 120)
        )
    params:
        gtf = "resources/Homo_sapiens.gtf",
        sjdbOverhang = "100",
        extra = ""
    log:
        "logs/star/index.log"
    wrapper:
        f"{swv}/bio/star/index"


rule star_mapping:
    input:
        unpack(sample_pair_w),
        index = "star/index/Homo_sapiens"
    output:
        temp("star/bam/{sample}/Aligned.sortedByCoord.out.bam")
    message:
        "Mapping {wildcards.sample} with STAR"
    threads:
        10
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(35840 + 10240 * attempt, 51200)
        ),
        time_min = (
            lambda wildcards, attempt: min(50 * attempt, 120)
        )
    params:
        index = "star/index/Homo_sapiens",
        extra = ("--outSAMtype BAM SortedByCoordinate "
                 "--outSAMattributes All "
                 "--outFilterType BySJout "
                 "--outFilterMismatchNmax 999 "
                 "--alignSJDBoverhangMin 1 "
                 "--outFilterMismatchNoverReadLmax 0.04 "
                 "--alignMatesGapMax 1000000 "
                 "--alignIntronMax 1000000 "
                 "--alignIntronMin 20 "
                 "--alignSJoverhangMin 8 "
                 "--outFilterMultimapNmax 20 "
                 "--twopassMode Basic")
    log:
        "logs/star/bam/{sample}.log"
    wrapper:
        f"{swv}/bio/star/align"


rule star_rename:
    input:
        "star/bam/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "star/bam/{sample}.bam"
    message:
        "Renaming {wildcards.sample} for further analyses"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(128 * attempt, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(5 * attempt, 15)
        )
    params:
        ""
    conda:
        "../envs/bash.yaml"
    log:
        "logs/rename/{sample}.log"
    shell:
        "cp -v {input} {output} > {log} 2>&1"


rule samtools_index:
    input:
        "star/bam/{sample}.bam"
    output:
        "star/bam/{sample}.bam.bai"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(512 * attempt, 1024)
        ),
        time_min = (
            lambda wildcards, attempt: min(15 * attempt, 45)
        )
    params:
        ""
    wrapper:
        f"{swv}/bio/samtools/index"


rule samtools_flagstat:
    input:
        "star/bam/{sample}.bam"
    output:
        temp("samtools/flagstat/{sample}.flagstat")
    message:
        "Gathering statistics on {wildcards.sample}'s mapping"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(512 * attempt, 1024)
        ),
        time_min = (
            lambda wildcards, attempt: min(15 * attempt, 45)
        )
    params:
        ""
    wrapper:
        f"{swv}/bio/samtools/flagstat"


rule download_fasta_cdna:
    output:
        temp("resources/Homo_sapiens_cdna.fasta")
    message:
        "Downloading Fasta from ensembl"
    params:
        species = "homo_sapiens",
        release = "98",
        datatype = "cdna",
        build = "GRCh38"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(128 * attempt, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(25 * attempt, 120)
        )
    log:
        "logs/ensembl_annotation/get_genome.log"
    wrapper:
        f"{swv}/bio/reference/ensembl-sequence"


rule correct_cdna_names:
    input:
        "resources/Homo_sapiens_cdna.fasta"
    output:
        temp("resources/Homo_sapiens_cdna_corrected.fasta")
    message:
        "Correcting CDNA names"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    conda:
        "../../envs/bash.yaml"
    log:
        "logs/correct_cdna_names.log"
    shell:
        "sed 's/\.[0-9]\+ / /g' {input} > {output} 2> {log}"


rule salmon_index:
    input:
        "resources/Homo_sapiens_cdna_corrected.fasta"
    output:
        directory("salmon/index/Homo_sapiens")
    message:
        "Indexing hsa with Salmon"
    threads:
        10
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(8192 + 1024 * attempt, 15360)
        ),
        time_min = (
            lambda wildcards, attempt: min(25 * attempt, 120)
        )
    params:
        extra = "--keepDuplicates --perfectHash"
    log:
        "logs/salmon/index.log"
    wrapper:
        f"{swv}/bio/salmon/index"


rule salmon_quant:
    input:
        unpack(salmon_sample_pair_w),
        index = "salmon/index/Homo_sapiens",
        gtf = ancient("resources/Homo_sapiens.gtf")
    output:
        quant = "salmon/quant/{sample}/quant.sf",
        genes = "salmon/quant/{sample}/quant.genes.sf"
    message:
        "Quantifying {wildcards.sample} with Salmon"
    threads:
        20
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(10240 + 1024 * attempt, 15360)
        ),
        time_min = (
            lambda wildcards, attempt: min(120 * attempt, 240)
        )
    params:
        libtype = "A",
        extra = "--numBootstraps 100 --validateMappings --gcBias --seqBias --geneMap resources/Homo_sapiens.gtf"
    log:
        "logs/salmon/quant/{sample}.log"
    wrapper:
        f"{swv}/bio/salmon/quant"


rule multiqc:
    input:
        expand(
            "salmon/quant/{sample}/quant.sf",
            sample=design.index.tolist()
        ),
        # expand(
        #     "samtools/flagstat/{sample}.flagstat",
        #     sample=design.index.tolist()
        # ),
        expand(
            "qc/fastqc/{sample}{ext}",
            sample=samplsdict.keys(),
            ext=[".html", "_fastqc.zip"]
        )
    output:
        report(
            "qc/multiqc.html",
            category="Quality Reports",
            caption="../../reports/MultiQC.report.rst"
        )
    message:
        "Gathering quality controls"
    params:
        ""
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(1024 + 1024 * attempt, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(30 * attempt, 120)
        )
    log:
        "logs/multiqc.log"
    wrapper:
        f"{swv}/bio/multiqc"

rule tx2gene:
    input:
        gtf = "resources/Homo_sapiens.gtf"
    output:
        tsv = temp("deseq2/Homo_sapiens.tsv")
    message:
        "Building T2G table"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 2048
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 119
        )
    params:
        gencode = True
    log:
        "logs/tx2gene.log"
    wrapper:
        f"{git}/tx_to_tgene/bio/tx_to_gene/gtf"


rule tr2gene:
    input:
        "deseq2/Homo_sapiens.tsv"
    output:
        temp("deseq2/Homo_sapiens.tximport.tsv")
    message:
        "Reducing the tx2gene table for tximport"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512, 1024)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 10, 20)
        )
    log:
        "logs/tx2gene.log"
    shell:
        "cut -f1,3 {input} > {output} 2> {log}"


rule tximport:
    input:
        quant = expand(
            "salmon/quant/{sample}/quant.sf",
            sample=design.index.tolist()
        ),
        tx_to_gene = "deseq2/Homo_sapiens.tximport.tsv"
    output:
        txi = "deseq/tximport.RDS"
    message:
        "Importing salmon data into R"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 20480
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 119
        )
    params:
        extra = "type = 'salmon'"
    log:
        "logs/tximport.log"
    wrapper:
        f"{swv}/bio/tximport"


rule pandas_merge:
    input:
        quant = expand(
            "salmon/quant/{sample}/quant.genes.sf",
            sample=design.index.tolist()
        ),
        tx2gene = "deseq2/Homo_sapiens.tsv"
    output:
        tsv = report(
            "salmon/aggregated/TPM.counts.tsv",
            caption="../../reports/aggregated.TPM.counts.rst",
            category="Counts"
        )
    message:
        "Aggregating salmon counts"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    params:
        drop_null = True,
        drop_na = True,
        gencode = True
    log:
        "logs/pandas_merge.log"
    wrapper:
        f"{git}/pandas-merge/bio/pandas/salmon"


rule box_count:
    input:
        counts = "salmon/aggregated/TPM.counts.tsv"
    output:
        png = report(
            "figures/box_counts.png",
            caption="../../reports/box.counts.rst",
            category="Figures"
        )
    message:
        "Plotting comparative boxplots of each sample's counts"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    log:
        "logs/box_count.log"
    wrapper:
        f"{git}/pandas-merge/bio/seaborn/box-counts"


rule pairwise_scatterplot:
    input:
        counts = "salmon/aggregated/TPM.counts.tsv"
    output:
        png = report(
            "figures/pairwise_scatterplot.png",
            caption="../../reports/pairwise_scatterplot.rst",
            category="Figures"
        )
    message:
        "Drawing pairwise scatterplot"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    log:
        "logs/pairwise_scatterplot.log"
    wrapper:
        f"{git}/pandas-merge/bio/seaborn/pairwise-scatterplot"


rule pca_plots:
    input:
        counts = "salmon/aggregated/TPM.counts.tsv"
    output:
        png = report(
            "figures/PCA/PCA_{factor}_PC1_PC2.png",
            caption="../../reports/pca.rst",
            category="Figures"
        )
    message:
        "Plotting PCA for factor {wildcards.factor}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 2048
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 40
        )
    params:
        axes = [1, 2],
        conditions = lambda wildcards : {
            k: v for k, v in zip(
                design.index.tolist(),
                design[wildcards.factor].tolist()
            )
        },
        prefix = (
            lambda wildcards: f"PCA_{wildcards.factor}"
        )
    log:
        "logs/pca_plots/{factor}.log"
    wrapper:
        f"{git}/pandas-merge/bio/seaborn/pca"


rule deseq2_dds:
    input:
        tximport = ancient("deseq/tximport.RDS"),
        coldata = design_path
    output:
        dds = temp("deseq2/RDS/DESeq2_dds_{factor}.RDS")
    message:
        "Building Deseq2 dataset for {wildcards.factor}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 40, 200)
        )
    params:
        design = (
            lambda wildcards: {
                "Liver": '~Liver',
                "Peritoneum": '~Peritoneum',

            }[wildcards.factor]
        )
    log:
        "logs/deseq2/dds/{factor}.log"
    wrapper:
        f"{git}/deseq2_dataset/bio/deseq2/DESeqDataSetFromTximport"


rule deseq2_esf:
    input:
        dds = "deseq2/RDS/DESeq2_dds_{factor}.RDS"
    output:
        dds = temp("deseq2/RDS/DEseq2_esf_{factor}.RDS")
    message:
        "Estimating size factors for {wildcards.factor}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 40, 200)
        )
    params:
        locfunc = "median"
    log:
        "logs/deseq2/esf/{factor}.log"
    wrapper:
        f"{git}/deseq2-estimateSizeFactors/bio/deseq2/estimateSizeFactors"


rule deseq2_disp:
    input:
        dds = "deseq2/RDS/DEseq2_esf_{factor}.RDS"
    output:
        disp = temp("deseq2/RDS/DESeq2_DispEst_{factor}.RDS")
    message:
        "Estimating sample dispersion for {wildcards.factor}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 200)
        )
    log:
        "logs/deseq2/disp/{factor}.log"
    params:
        fittype = "parametric"
    wrapper:
        f"{git}/deseq2-disp.R/bio/deseq2/estimateDispersions"


rule deseq2_vst:
    input:
        dds = "deseq2/RDS/DESeq2_DispEst_{factor}.RDS"
    output:
        rds = temp("deseq2/RDS/DESeq2_vst_{factor}.RDS"),
        tsv = "deseq2/TSV/DEseq2_VST_{factor}.tsv"
    message:
        "Computing variance stabilizing transform for {wildcards.factor}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 8192
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 40
        )
    log:
        "logs/deseq2/vst/{factor}.log"
    wrapper:
        f"{git}/deseq2-vst/bio/deseq2/vst"

rule deseq2_waldtest:
    input:
        dds = "deseq2/RDS/DESeq2_DispEst_{factor}.RDS"
    output:
        rds = "deseq2/RDS/DESeq2_WaldTest_{factor}.RDS",
        tsv = report(
            directory("deseq2/TSV/DESeq2_{factor}"),
            caption="../../reports/deseq2.complete.table.rst",
            category="Differential Expression"
        )
    message:
        "Performing DESeq2 wald test on {wildcards.factor}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 5125
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 30
        )
    log:
        "logs/deseq2_waldtest/{factor}.log"
    wrapper:
        f"{git}/deseq2-waldtest/bio/deseq2/nbinomWaldTest"


## TODO: Hierarchical sample clustering
## TODO: Factor correlation
## TODO: Immune decov
## TODO: Clustering on Immune decov
