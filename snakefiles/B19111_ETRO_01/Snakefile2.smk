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
container: "docker://continuumio/miniconda3:5.0.1"

design = pandas.read_csv(
    "../design.tsv",
    header=0,
    index=0
)

rule all:
    input:
        expand("bcftools/norm/{sample}.vcf.gz", sample=design.index)
    message:
        "Finishing pipeline"


rule rename:
    input:
        "/mnt/beegfs/scratch/bioinfo_core/B19111_ETRO_01/annovar/raw_data/{sample}.vcf"
    output:
        "raw_data/{sample}.vcf.gz"
    message:
        "Rename gzipped files"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 5, 200)
        )
    group:
        "prepare"
    conda:
        "../../envs/bash.yaml"
    params:
        "-s"
    log:
        "logs/link/{sample}"
    shell:
        "ln {params} {input} {output} 2> {log}"


rule tabix_index:
    input:
        "raw_data/{sample}.vcf.gz"
    output:
        "raw_data/{sample}.vcf.tbi"
    message:
        "Indexing {wildcards.sample}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    group:
        "prepare"
    params:
        "-p vcf"
    log:
        "logs/tabix/{sample}.log"
    wrapper:
        f"{git}/bio/tabix"


rule norm_vcf:
    input:
        "raw_data/{sample}.vcf.gz"
    output:
        "bcftools/norm/{sample}.vcf.gz"
    message:
        "Normalizing {wildcars.sample}"
    threads:
        2
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    params:
        "--multiallelics '-' --output-type z --threads 1"
    wrapper:
        f"{git}/bio/bcftools/norm"
