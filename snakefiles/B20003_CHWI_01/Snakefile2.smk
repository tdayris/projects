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

nb = [f"000{n}" for n in range(0, 5, 1)]
sample = ["S10", "S12", "S14", "S15"]
samples = ["S2", "S3", "S11", "S13"] + sample

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


wildcard_constraints:
    sample = "|".join(sample),
    nb = "|".join(nb),
    samples = "|".join(samples)


rule all:
    input:
        "qc/Variant_statistics.html",
        "JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.tsv",
        expand("tables/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.tsv", sample=sample)
    message:
        "Finishing pipeline"


rule compression:
    input:
        vcf = "bcftools/call/{samples}.vcf"
    output:
        vcf = "bcftools/compress/{samples}.vcf.gz",
        index = "bcftools/compress/{samples}.vcf.tbi"
    message:
        "Comprssing and indexing {wildcards.samples}"
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
    params:
        bgzip = "-c",
        tabix = "-p vcf"
    log:
        compress = "logs/bgzip/bcftools_call_{samples}.log",
        tabix = "logs/tabix/bcftools_call_{samples}.log"
    shell:
        " bgzip {params.bgzip} {input.vcf} "
        " > {output.vcf} 2> {log.compress} "
        " && "
        " tabix {params.tabix} {output.vcf} "
        " > {log.tabix} 2>&1 "


rule define_base:
    input:
        S2 = "bcftools/compress/S2.vcf.gz",
        S2_index = "bcftools/compress/S2.vcf.gz.tbi",
        S3 = "bcftools/compress/S3.vcf.gz",
        S3_index = "bcftools/compress/S3.vcf.gz.tbi"
    output:
        baseline = "bcftools/isec/JAK2.vcf.gz"
    message:
        "Building baseline with S2 and S3"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 4096, 20480)
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 15
        )
    log:
        "logs/bcftools/isec/baseline.log"
    conda:
        "../../envs/biotools.yaml"
    params:
        isec = (
            " --collapse none "
            " --exclude '(INFO/DP < 40)' "
            " --output-type z "
            " --threads 1 "
            " -n '-1' "
        )
    shell:
        " bcftools isec "
        " {params.isec} "
        " {input.S2} "
        " {input.S3} "
        " > {output.baseline} "
        " 2> {log} "


rule snpeff_baseline:
    input:
        "bcftools/isec/JAK2.vcf.gz"
    output:
        calls="snpeff/JAK2/JAK2.vcf",
        stats="snpeff/JAK2/JAK2.html",
        csvstats="snpeff/JAK2/JAK2.csv"
    message:
        "Annotating JAK2 with snpeff"
    threads:
        2
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


rule define_no_prolif:
    input:
        S2 = "bcftools/compress/S2.vcf.gz",
        S2_index = "bcftools/compress/S2.vcf.gz.tbi",
        S3 = "bcftools/compress/S3.vcf.gz",
        S3_index = "bcftools/compress/S3.vcf.gz.tbi",
        S10 = "bcftools/compress/S10.vcf.gz",
        S10_index = "bcftools/compress/S10.vcf.gz.tbi",
        S11 = "bcftools/compress/S11.vcf.gz",
        S11_index = "bcftools/compress/S11.vcf.gz.tbi"
    output:
        no_prolif = "bcftools/isec/JAK2_vs_JAK2_SRSF2_S10_S11.vcf.gz"
    message:
        "Building no-proliferation mask with S2 and S3 in one hand, "
        "and S10 and S11 in the other hand"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 4096, 20480)
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 15
        )
    log:
        "logs/bcftools/isec/baseline.log"
    conda:
        "../../envs/biotools.yaml"
    params:
        isec = (
            " --collapse none "
            " --exclude '(INFO/DP < 40)' "
            " --output-type z "
            " --threads 1 "
            " -n '~0011' "
        )
    shell:
        " bcftools isec "
        " {params.isec} "
        " {input.S2} "
        " {input.S3} "
        " {input.S10} "
        " {input.S11} "
        " > {output.no_prolif} "
        " 2> {log} "



rule snpeff_no_prolif:
    input:
        "bcftools/isec/JAK2_vs_JAK2_SRSF2_S10_S11.vcf.gz"
    output:
        calls="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_S10_S11.vcf",
        stats="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_S10_S11.html",
        csvstats="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_S10_S11.csv"
    threads:
        2
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



rule bcftools_isec:
    input:
        S2 = "bcftools/compress/S2.vcf.gz",
        S2_index = "bcftools/compress/S2.vcf.gz.tbi",
        S3 = "bcftools/compress/S3.vcf.gz",
        S3_index = "bcftools/compress/S3.vcf.gz.tbi",
        S10 = "bcftools/compress/S10.vcf.gz",
        S10_index = "bcftools/compress/S10.vcf.gz.tbi",
        S11 = "bcftools/compress/S11.vcf.gz",
        S11_index = "bcftools/compress/S11.vcf.gz.tbi",
        SXX = "bcftools/compress/{sample}.vcf.gz",
        SXX_index = "bcftools/compress/{sample}.vcf.gz.tbi"
    output:
        comparison = expand(
            "bcftools/isec/JAK2_vs_JAK2_SRSF2_{sample}/{nb}.vcf.gz",
            nb = nb,
            allow_missing = True
        )
    message:
        "Building no-proliferation mask with S2 and S3 in one hand, "
        "and S10 and S11 in the other hand"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 4096, 20480)
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 15
        )
    log:
        "logs/bcftools/isec/{sample}.log"
    conda:
        "../../envs/biotools.yaml"
    params:
        isec = lambda wildcards: f" --collapse none --exclude '(INFO/DP < 40)' --output-type z --threads 1 --prefix bcftools/isec/JAK2_vs_JAK2_SRSF2_{wildcards.sample}/ -n '~00001'"
    shell:
        " bcftools isec "
        " {params.isec} "
        " {input.S2} "
        " {input.S3} "
        " {input.S10} "
        " {input.S11} "
        " {input.SXX} "
        " > {log} "
        " 2>&1 "


rule bcftools_rename:
    input:
        vcf = "bcftools/isec/JAK2_vs_JAK2_SRSF2_{sample}/0004.vcf.gz"
    output:
        vcf = "bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz"
    message:
        "Renaming vcf for {wildcards.sample}"
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
    params:
        "--verbose"
    log:
        "logs/bcftools/rename/{sample}.log"
    shell:
        " cp {params} {input.vcf} {output.vcf} > {log} 2>&1 "


rule snpeff_compare:
    input:
        "bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz"
    output:
        calls="snpeff/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf",
        stats="snpeff/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.html",
        csvstats="snpeff/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.csv"
    message:
        "Annotating {wildcards.sample} with snpeff"
    threads:
        2
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
        "logs/snpeff/variants_present_only_in_{sample}.log"
    wrapper:
        f"{git}/bio/snpeff"


rule multiqc:
    input:
        expand(
            "snpeff/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.html",
            sample=sample
        ),
        expand(
            "snpeff/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.csv",
            sample=sample
        ),
        "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_S10_S11.html",
        "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_S10_S11.csv",
        "snpeff/JAK2/JAK2.html",
        "snpeff/JAK2/JAK2.csv",
        "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.html",
        "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.csv"
    output:
        report(
            "qc/Variant_statistics.html",
            caption="../../reports/MultiQC.report.rst",
            category="Quality Reports"
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
        "logs/multiqc.log"
    wrapper:
        f"{git}/bio/multiqc"


rule vcf_to_tsv_snpsift_prolif:
    input:
        vcf = "bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz"
    output:
        tsv = report(
            "tables/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.tsv",
            caption="../../reports/tsv.rst",
            category="Results"
        )
    message:
        "Converting VCF to TSV for {wildcards.sample} only"
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
    params:
        fields = " ".join(TSV_fields),
        extra = "-e '.' -s ','"
    log:
        "logs/snpsift/table/{sample}.log"
    shell:
        "SnpSift extractFields "
            " {params.extra} "
            " {input.vcf} "
        " {params.fields} "
        " > {output.tsv} "
        " 2> {log}"


rule prolif_common:
    input:
        expand("bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz", sample=sample)
    output:
        "bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.vcf.gz"
    message:
        "Finding common variants among samples with a proliferation relapse"
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
    params:
        isec = lambda wildcards: (
            " --collapse none "
            " --exclude '(INFO/DP < 40)' "
            " --output-type z "
            " --threads 1 "
            " -n '1111' "
        )
    log:
        "logs/bcftools/isec/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.log"
    shell:
        " bcftools isec "
        " {params.isec} "
        " {input} "
        " > {output} "
        " 2> {log} "


rule snpeff_common:
    input:
        "bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.vcf.gz"
    output:
        calls="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.vcf",
        stats="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.html",
        csvstats="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.csv"
    message:
        "Annotating common variants among proliferation with snpeff"
    threads:
        2
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
        "logs/snpeff/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.log"
    wrapper:
        f"{git}/bio/snpeff"


rule vcf_to_tsv_snpsift_common:
    input:
        vcf = "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.vcf"
    output:
        tsv = report(
            "JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.tsv",
            caption="../../reports/tsv.rst",
            category="Results"
        )
    message:
        "Converting VCF to TSV for JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15 only"
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
    params:
        fields = " ".join(TSV_fields),
        extra = "-e '.' -s ','"
    log:
        "logs/snpsift/table/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.log"
    shell:
        "SnpSift extractFields "
            " {params.extra} "
            " {input.vcf} "
        " {params.fields} "
        " > {output.tsv} "
        " 2> {log}"
