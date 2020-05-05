
import sys
import pandas
import os
import os.path

from snakemake.utils import min_version, makedirs
from pathlib import Path
from typing import Any, Dict

if sys.version_info < (3, 8):
    raise SystemError("Please use Python 3.8 or later")

min_version('5.16.0')
wrapper_version = '0.51.0'
git = "https://raw.githubusercontent.com/tdayris/snakemake-wrappers/Unofficial"
singularity: "docker://continuumio/miniconda3:5.0.1"

QUANT = config["quant_path"]
salmon_dirs = config["rna-count-salmon"]
workdir: config["workdir"]

SAMPLES = [
    os.path.basename(i)
    for i in "pseudo_mapping/quant/"
    if os.path.isdir(i)
]

sample_constraint = "|".join(SAMPLES)


rule all:
    input:
        quantiseq = "quanTIseq/results"
    message:
        "Finishing pipeline "



rule filter_quant:
    input:
        QUANT
    output:
        "immunedeconv/TPM.tsv"
    message:
        "Filtering the TPM table to fit immunedeconv requirements"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512, 1024)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 5, 15)
        )
    conda:
        "../../envs/bash.yaml"
    log:
        "logs/filter_quant.log"
    shell:
        "(paste <(head -n1 {input} | cut -f 35) <(head -n1 {input} | cut -f 2-34) && "
        "paste <(cut -f 35 {input}) <(cut -f 2-34 {input}) | "
        "sed '1d' "
        "| sort "
        "| uniq -w 7 -u"
        "| sed 's/Hugo_ID/GENE/g' "
        "| grep -vP \"^\\t\")"
        "> {output} 2> {log}"

rule quanTIseq_pipeline:
    input:
        "immunedeconv/TPM.tsv"
    output:
        directory("quanTIseq/results")
    message:
        "Performing quanTIseq estimation"
    threads:
        10
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    conda:
        "../../envs/immunedeconv.yaml"
    params:
        path = "/mnt/beegfs/userdata/t_dayris/projects",
        extra = (
            "--pipelinestart=decon "
            "--tumor=TRUE "
            "--mRNAscale=TRUE "
        )
    log:
        "logs/quanTIseq_pipeline.log"
    shell:
        "{params.path}/scripts/quanTIseq_pipeline.sh "
        "--inputfile={input} "
        "--outputdir={output} "
        "--threads={threads} "
        "{params.extra} "
        "> {log} 2>&1"


rule quanTIseq_deconv:
    input:
        quant = "immunedeconv/TPM.tsv"
    output:
        rds = "immunedeconv/quantiseq.RDS"
    message:
        "Performing quantiseq"
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
        "envs/immunedeconv.yaml"
    log:
        "logs/quanTIseq.log"
    script:
        "../../scripts/immunedeconv_quantiseq.R"
