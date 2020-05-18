import sys
import pandas
import os.path

from snakemake.utils import min_version, makedirs
from pathlib import Path
from typing import Any, Dict

if sys.version_info < (3, 8):
    raise SystemError("Please use Python 3.8 or later")

min_version('5.16.0')
git = "https://raw.githubusercontent.com/tdayris/snakemake-wrappers/Unofficial"
singularity: "docker://continuumio/miniconda3:5.0.1"

samples = {
    "S2": "JAK2_S2",
    "S3": "JAK2_S3",
    "S10": "JAK2_SRSF2_S10",
    "S11": "JAK2_SRSF2_S11",
    "S12": "JAK2_SRSF2_S12",
    "S13": "JAK2_SRSF2_S13",
    "S14": "JAK2_SRSF2_S14",
    "S15": "JAK2_SRSF2_S15",
}
rev_sample = {v: k for k, v in samples.items()}
numbers_list = [f"000{i}" for i in range(1, 4, 1)]


TSV_fields = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "GT",
    "ANN[*].GENE",
    "ANN[*].EFFECT",
    "ANN[*].IMPACT",
    "ANN[*].BIOTYPE",
    "ANN[*].HGVS_C",
    "ANN[*].HGVS_P",
    "ANN[*].CDNA_POS",
    "DP",
    "ANN[*].GENEID",
]

report: "general.rst"
wildcard_constraints:
    sample_short_id = "|".join(samples.keys()),
    sample = "|".join(samples.values()),
    numbers = "|".join(numbers_list)

design_path = Path("design_calling.csv")
if not design_path.exists():
    raise FileNotFoundError(f"Could not find {str(design_path)}")

design = pandas.read_csv(
    design_path,
    sep="\t",
    header=0,
    index_col=0,
    dtype=str
)
design_dict = design.to_dict()
FASTA = "resources/Mus_musculus.GRCm38.dna.toplevel.fa"

rule all:
    input:
        "qc/Complete_report.html",
        "qc/Baseline_report.html"
    message:
        "Finishing pipeline"


rule bcftools_call:
    input:
        ref = FASTA,
        samples = lambda wildcards: design_dict["Bam_File"][wildcards.sample],
        indexes = lambda wildcards: design_dict["Bam_Index"][wildcards.sample]
    output:
        "bcftools/call/{sample}.vcf.gz"
    message:
        "Calling variants on {wildcards.sample} with BCFTools"
    threads:
        3
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 20480)
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 115
        )
    params:
        mpileup = "--adjust-MQ 50 --min-MQ 10 --min-BQ 15",
        call = "--prior 0.001 --variants-only --output-type z --threads 1"
    log:
        "logs/bcftools/call/{sample}.log"
    wrapper:
        f"{git}/bio/bcftools/call"


rule tabix_call:
    input:
        "bcftools/call/{sample}.vcf.gz"
    output:
        "bcftools/call/{sample}.vcf.gz.tbi"
    message:
        "Indexing snpeff result for {wildcards.sample}"
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
        '-p vcf'
    log:
        "logs/tabix/snpeff_{sample}.log"
    wrapper:
        f"{git}/bio/tabix"


rule define_control_baseline:
    input:
        S2 = "bcftools/call/JAK2_S2.vcf.gz",
        S2_index = "bcftools/call/JAK2_S2.vcf.gz.tbi",
        S3 = "bcftools/call/JAK2_S3.vcf.gz",
        S3_index = "bcftools/call/JAK2_S3.vcf.gz.tbi"
    output:
        merge = "comparisons/baseline.vcf"
    message:
        "Building control baseline"
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
        "../../envs/biotools.yaml"
    log:
        "logs/bcftools/merge/define_control_baseline.log"
    params:
        "--merge none --output-type v"
    shell:
        " bcftools merge "
        " {params} "
        " {input.S2} "
        " {input.S3} "
        " --output {output.merge} "
        " > {log} 2>&1 "


rule define_no_prolif_baseline:
    input:
        S11 = "bcftools/call/JAK2_SRSF2_S11.vcf.gz",
        S11_index = "bcftools/call/JAK2_SRSF2_S11.vcf.gz.tbi",
        S13 = "bcftools/call/JAK2_SRSF2_S13.vcf.gz",
        S13_index = "bcftools/call/JAK2_SRSF2_S13.vcf.gz.tbi"
    output:
        merge = "comparisons/no_prolif_baseline.vcf"
    message:
        "Building baseline for mutant with no proliferation"
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
        "../../envs/biotools.yaml"
    log:
        "logs/bcftools/merge/define_no_prolif_baseline.log"
    params:
        "--merge none --output-type v"
    shell:
        " bcftools merge "
        " {params} "
        " {input.S11} "
        " {input.S13} "
        " --output {output.merge} "
        " > {log} 2>&1 "


rule overlap_ctrl_and_no_prolif:
    input:
        baseline = "comparisons/baseline.vcf",
        no_prolif = "comparisons/no_prolif_baseline.vcf"
    output:
        isec = expand(
            "comparison/overlap_ctrl_and_no_prolif/{subsets}.vcf",
            subsets = [
                "JAK2_only",
                "JAK2_SRSF2_only",
                "JAK2_also_in_JAK2_SRSF2",
                "JAK2_SRSF2_also_in_JAK2"
            ]
        )
    message:
        "All overlapping sets for control and non-proliferation samples"
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
        "../../envs/biotools.yaml"
    log:
        isec = "logs/overlap_ctrl_and_no_prolif.log",
        mv_jak2_only = "logs/mv_jak2_only.log",
        mv_jak2_shared = "logs/mv_jak2_shared.log",
        mv_jak2_srsf2_only = "logs/mv_jak2_srsf2_only.log",
        mv_jak2_srsf2_shared = "logs/mv_jak2_srsf2_shared.log"
    params:
        isec = "--collapse none --exclude '(INFO/DP < 40)' --output-type v",
        mv = "--verbose --force",
        prefix = "comparison/overlap_ctrl_and_no_prolif"
    shell:
        " bcftools isec "
        " {params.isec} "
        " {input.baseline} "
        " {input.no_prolif} "
        " -p {params.prefix} "
        " > {log.isec} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0000.vcf "
        " {params.prefix}/JAK2_only.vcf "
        " > {log.mv_jak2_only} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0001.vcf "
        " {params.prefix}/JAK2_also_in_JAK2_SRSF2.vcf "
        " > {log.mv_jak2_shared} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0002.vcf "
        " {params.prefix}/JAK2_SRSF2_only.vcf "
        " > {log.mv_jak2_srsf2_only} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0003.vcf "
        " {params.prefix}/JAK2_SRSF2_also_in_JAK2.vcf "
        " > {log.mv_jak2_srsf2_shared} 2>&1 "


rule define_ctrl_and_no_prolif_baseline:
    input:
        S11 = "bcftools/call/JAK2_SRSF2_S11.vcf.gz",
        S11_index = "bcftools/call/JAK2_SRSF2_S11.vcf.gz.tbi",
        S13 = "bcftools/call/JAK2_SRSF2_S13.vcf.gz",
        S13_index = "bcftools/call/JAK2_SRSF2_S13.vcf.gz.tbi",
        S2 = "bcftools/call/JAK2_S2.vcf.gz",
        S2_index = "bcftools/call/JAK2_S2.vcf.gz.tbi",
        S3 = "bcftools/call/JAK2_S3.vcf.gz",
        S3_index = "bcftools/call/JAK2_S3.vcf.gz.tbi"
    output:
        merge = "comparisons/ctrl_and_no_prolif_baseline.vcf"
    message:
        "Building baseline for control and mutants without proliferation"
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
        "../../envs/biotools.yaml"
    log:
        "logs/bcftools/merge/ctrl_and_no_prolif_baseline.log"
    params:
        "--merge none --output-type v"
    shell:
        " bcftools merge "
        " {params} "
        " {input.S2} "
        " {input.S3} "
        " {input.S11} "
        " {input.S13} "
        " --output {output.merge} "
        " > {log} 2>&1 "


rule compress_and_index_basline:
    input:
        vcf = "comparisons/ctrl_and_no_prolif_baseline.vcf"
    output:
        vcf = "comparisons/ctrl_and_no_prolif_baseline.vcf.gz",
        tbi = "comparisons/ctrl_and_no_prolif_baseline.vcf.gz.tbi"
    message:
        "Indexing and compressing baseline"
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
        "../../envs/biotools.yaml"
    log:
        bgzip = "logs/compress_and_index_basline/bgzip.log",
        tabix = "logs/compress_and_index_basline/tabix.log",
    shell:
        "bgzip -c {input.vcf} > {output.vcf} 2> {log.bgzip} && "
        " tabix -p vcf {output.vcf} > {log.tabix} 2>&1 "


rule overlap_ctrl_and_prolif:
    input:
        baseline = "comparisons/ctrl_and_no_prolif_baseline.vcf.gz",
        baseline_idx = "comparisons/ctrl_and_no_prolif_baseline.vcf.gz.tbi",
        no_prolif = "bcftools/call/{sample}.vcf.gz",
        no_prolif_index = "bcftools/call/{sample}.vcf.gz.tbi",
    output:
        baseline_only = "comparison/prolif_{sample}/baseline_S2_S3_S11_S13_not_in_{sample}.vcf",
        baseline_shared = "comparison/prolif_{sample}/baseline_S2_S3_S11_S13_shared_with_{sample}.vcf",
        sample_only = "comparison/prolif_{sample}/{sample}_only.vcf",
        sample_shared = "comparison/prolif_{sample}/{sample}_shared.vcf",
    message:
        "All overlapping sets for control and "
        "non-proliferation vs {wildcards.sample}"
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
        "../../envs/biotools.yaml"
    log:
        isec = "logs/overlap_ctrl_and_no_prolif_{sample}.log",
        mv_baseline_only = "logs/mv_baseline_only_{sample}.log",
        mv_baseline_shared = "logs/mv_baseline_shared_{sample}.log",
        mv_sample_only = "logs/mv_{sample}_only.log",
        mv_sample_shared = "logs/mv_{sample}_shared.log"
    params:
        isec = "--collapse none --exclude '(INFO/DP < 40)' --output-type v",
        mv = "--verbose --force",
        prefix = lambda wildcards: f"comparison/prolif_{wildcards.sample}"
    shell:
        " bcftools isec "
        " {params.isec} "
        " {input.baseline} "
        " {input.no_prolif} "
        " -p {params.prefix} "
        " > {log.isec} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0000.vcf "
        " {output.baseline_only} "
        " > {log.mv_baseline_only} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0001.vcf "
        " {output.baseline_shared} "
        " > {log.mv_baseline_shared} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0002.vcf "
        " {output.sample_only} "
        " > {log.mv_sample_only} 2>&1 && "
        " mv {params.mv} "
        " {params.prefix}/0003.vcf "
        " {output.sample_shared} "
        " > {log.mv_sample_shared} 2>&1 "


rule baseline_annotate:
    input:
        "comparisons/ctrl_and_no_prolif_baseline.vcf"
    output:
        "comparisons/snpeff/ctrl_and_no_prolif_baseline.vcf"
    message:
        "Annotating JAK2 and JAK2-SRSF2-no-prolif with snpeff"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 4096, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 30, 200)
        )
    params:
        reference = "GRCm38.86",
        extra = "-Xmx4g"
    log:
        "logs/snpeff/JAK2.log"
    wrapper:
        f"{git}/bio/snpeff"


rule annotate_prolif:
    input:
        "comparison/prolif_{sample}/{sample}_{status}.vcf"
    output:
        calls = "comparison/snpeff/call/prolif_{sample}/{sample}_{status}.vcf",
        stats = "comparison/snpeff/stats/prolif_{sample}/{sample}_{status}.csv",
        csvstats = "comparison/snpeff/report/prolif_{sample}/{sample}_{status}.tml"
    message:
        "Annotating {wildcards.sample} ({wildcards.status}) with snpeff"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 4096, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 30, 200)
        )
    params:
        reference = "GRCm38.86",
        extra = "-Xmx4g"
    log:
        "logs/snpeff/annotate_prolif_{sample}_{status}.log"
    wrapper:
        f"{git}/bio/snpeff"


rule annotate_prolif_baseline:
    input:
        "comparison/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.vcf"
    output:
        calls = "comparison/snpeff/call/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.vcf",
        stats = "comparison/snpeff/stats/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.csv",
        csvstats = "comparison/snpeff/report/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.html"
    message:
        "Annotating not-in-{wildcards.sample} ({wildcards.status}) with snpeff"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 4096, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 30, 200)
        )
    params:
        reference = "GRCm38.86",
        extra = "-Xmx4g"
    log:
        "logs/snpeff/annotate_prolif_baseline_{status}_{sample}.log"
    wrapper:
        f"{git}/bio/snpeff"


rule prolif_baseline_report:
    input:
        expand(
            "comparison/snpeff/stats/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.stats",
            sample=["JAK2_SRSF2_S10", "JAK2_SRSF2_S12", "JAK2_SRSF2_S14", "JAK2_SRSF2_S15"],
            status=["not_in", "shared_with"]
        ),
        expand(
            "comparison/snpeff/report/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.html",
            sample=["JAK2_SRSF2_S10", "JAK2_SRSF2_S12", "JAK2_SRSF2_S14", "JAK2_SRSF2_S15"],
            status=["not_in", "shared_with"]
        )
    output:
        report(
            "qc/Baseline_report.html",
            caption="reports/baseline_report.rst",
            category="Histograms"
        )
    message:
        "Compiling variant metrics"
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
        "logs/multiqc/baseline/baseline.log"
    wrapper:
        f"{git}/bio/multiqc"


rule complete_report:
    input:
        expand(
            "comparison/snpeff/stats/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.csv",
            sample=["JAK2_SRSF2_S10", "JAK2_SRSF2_S12", "JAK2_SRSF2_S14", "JAK2_SRSF2_S15"],
            status=["not_in", "shared_with"]
        ),
        expand(
            "comparison/snpeff/report/prolif_{sample}/baseline_S2_S3_S11_S13_{status}_{sample}.html",
            sample=["JAK2_SRSF2_S10", "JAK2_SRSF2_S12", "JAK2_SRSF2_S14", "JAK2_SRSF2_S15"],
            status=["not_in", "shared_with"]
        ),
        expand(
            "comparison/snpeff/report/prolif_{sample}/{sample}_{status}.html",
            sample=["JAK2_SRSF2_S10", "JAK2_SRSF2_S12", "JAK2_SRSF2_S14", "JAK2_SRSF2_S15"],
            status=["only", "shared"]
        ),
        expand(
            "comparison/snpeff/stats/prolif_{sample}/{sample}_{status}.csv",
            sample=["JAK2_SRSF2_S10", "JAK2_SRSF2_S12", "JAK2_SRSF2_S14", "JAK2_SRSF2_S15"],
            status=["only", "shared"]
        )
    output:
        report(
            "qc/Complete_report.html",
            caption="reports/complete_report.rst",
            category="Histograms"
        )
    message:
        "Compiling variant metrics"
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
        "logs/multiqc/complete.log"
    wrapper:
        f"{git}/bio/multiqc"
