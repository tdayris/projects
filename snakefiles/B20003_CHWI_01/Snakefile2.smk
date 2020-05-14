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


wildcard_constraints:
    sample = "|".join(sample),
    nb = "|".join(nb),
    samples = "|".join(samples)


rule all:
    input:
        baseline = "bcftool/isec/baseline.vcf.gz",
        no_prolif = "bcftool/isec/JAK2_vs_JAK2_SRSF2_S10_S11.vcf.gz",
        comparison = expand(
            "bcftool/isec/JAK2_vs_JAK2_SRSF2_{sample}/{nb}.vcf.gz",
            nb = nb,
            sample = sample
        )
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
        baseline = "bcftool/isec/baseline.vcf.gz"
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
        no_prolif = "bcftool/isec/JAK2_vs_JAK2_SRSF2_S10_S11.vcf.gz"
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
            "bcftool/isec/JAK2_vs_JAK2_SRSF2_{sample}/{nb}.vcf.gz",
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
        isec = (
            " --collapse none "
            " --exclude '(INFO/DP < 40)' "
            " --output-type z "
            " --threads 1 "
            " --prefix bcftool/isec/JAK2_vs_JAK2_SRSF2/ "
	    " -n '~00001'"
        )
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
