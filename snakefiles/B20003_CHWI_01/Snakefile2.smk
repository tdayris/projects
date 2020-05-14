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

nb = [f"000{n}" for n in range(0, 4, 1)]
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

report: "general.rst"

wildcard_constraints:
    sample = "|".join(sample),
    nb = "|".join(nb),
    samples = "|".join(samples)


rule all:
    input:
        "qc/Variant_statistics.html",
        "JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.tsv",
        expand("tables/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.tsv", sample=sample),
        expand("qc/Common_Variants_JAK2_SRSF2_with_proliferation_by_{jak_srsf}.html", jak_srsf = sample)
    message:
        "Finishing pipeline"


rule compression:
    input:
        vcf = "bcftools/call/{samples}.vcf"
    output:
        vcf = "bcftools/compress/{samples}.vcf.gz",
        index = "bcftools/compress/{samples}.vcf.gz.tbi"
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


rule multiqc_only:
    input:
        expand(
            "snpeff/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.html",
            sample=sample
        ),
        expand(
            "snpeff/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.csv",
            sample=sample
        )
    output:
        report(
            "qc/Unique_variants_in_{sample}.html",
            caption="../../reports/MultiQC.report.rst",
            category="Quality Reports"
        )
    message:
        "Compiling variant metrics for {wildcards.sample}"
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
        "logs/multiqc/JAK2_SRSF2/{sample}_unique.log"
    wrapper:
        f"{git}/bio/multiqc"



rule multiqc_common:
    input:
        expand(
            "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.html",
            sample=sample
        ),
        expand(
            "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.csv",
            sample=sample
        )
    output:
        report(
            "qc/Common_Variants_JAK2_SRSF2_with_proliferation_by_{sample}.html",
            caption="../../reports/MultiQC.report.rst",
            category="Quality Reports"
        )
    message:
        "Compiling variant metrics for {wildcards.sample}"
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
        "logs/multiqci/JAK2_SRSF2/{sample}_common.log"
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


rule tabix_sample_only:
    input:
        "bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz"
    output:
        "bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz.tbi"
    message:
        "Indexing {wildcards.sample}-only variants"
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
        "logs/tabix_sample_only/{sample}.log"
    wrapper:
        f"{git}/bio/tabix"


rule prolif_common:
    input:
        vcf = expand(
            "bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz",
            sample=sample
        ),
        idx = expand(
            "bcftools/isec/JAK2_vs_JAK2_SRSF2/variants_present_only_in_{sample}.vcf.gz.tbi",
            sample=sample
        )
    output:
        expand(
            "bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{nb}.vcf.gz",
            nb = nb
        )
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
            " -n '=3' "
            " -p bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/"
        )
    log:
        "logs/bcftools/isec/JAK2_vs_JAK2_SRSF2_common_S10_S12_S14_S15.log"
    shell:
        " bcftools isec "
        " {params.isec} "
        " {input.vcf} "
        " > {log} "
        " 2>&1 "


rule rename_commons:
    input:
        i1 = "bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/0000.vcf.gz",
        i2 = "bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/0001.vcf.gz",
        i3 = "bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/0002.vcf.gz",
        i4 = "bcftools/isec/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/0003.vcf.gz"
    output:
        o1 = "bcftools/renamed/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/S10.vcf.gz",
        o2 = "bcftools/renamed/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/S12.vcf.gz",
        o3 = "bcftools/renamed/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/S14.vcf.gz",
        o4 = "bcftools/renamed/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/S15.vcf.gz"
    message:
        "Renaming common searched within JAK2 SRSF2"
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
        l1 = "logs/rename/commons/1.log",
        l2 = "logs/rename/commons/2.log",
        l3 = "logs/rename/commons/3.log",
        l4 = "logs/rename/commons/4.log"
    shell:
        " cp {params} {input.i1} {output.o1} > {log} 2>&1 && "
        " cp {params} {input.i2} {output.o2} > {log} 2>&1 && "
        " cp {params} {input.i3} {output.o3} > {log} 2>&1 && "
        " cp {params} {input.i4} {output.o4} > {log} 2>&1 "


rule snpeff_common:
    input:
        "bcftools/renamed/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.vcf.gz"
    output:
        calls="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.vcf",
        stats="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.html",
        csvstats="snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.csv"
    message:
        "Annotating common variants among proliferation for {wildcards.sample}"
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
        "logs/snpeff/JAK2_vs_JAK2_SRSF2_common/{sample}.log"
    wrapper:
        f"{git}/bio/snpeff"


rule vcf_to_tsv_snpsift_common:
    input:
        vcf = "snpeff/JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.vcf"
    output:
        tsv = report(
            "JAK2_vs_JAK2_SRSF2/JAK2_vs_JAK2_SRSF2_common/{sample}.tsv",
            caption="../../reports/tsv.rst",
            category="Results"
        )
    message:
        "Converting VCF to TSV for JAK2_vs_JAK2_SRSF2 {wildcards.sample}"
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
        "logs/snpsift/table/JAK2_vs_JAK2_SRSF2_common/{sample}.log"
    shell:
        "SnpSift extractFields "
        " {params.extra} "
        " {input.vcf} "
        " {params.fields} "
        " > {output.tsv} "
        " 2> {log}"
