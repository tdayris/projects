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
containers: "docker://continuumio/miniconda3:5.0.1"

bam_path = Path("bams")
fasta_path = Path("../wes-mapping/")
gnomad_path = Path("resources/")
bam = {
    f.name[:-13]: f
    for f in bam_path.iterdir()
    if f.name.endswith("bam")
}

bai = {
    f.name[:-13]: f
    for f in bam_path.iterdir()
    if f.name.endswith("bai")
}

print(bam)
print(bai)

wildcard_constraints:
    samples = "|".join(bams.keys())
    bai = "|".join(bai.keys())


rule all:
    input:
        # TODO: Fill this with final files
    message:
        "Finishing pipeline"


rule mutect2:
    input:
        fasta = fasta_path.absolute(),
        map = lambda wildcards: bam[wildcards.sample],
        indexes = lambda wildcards: bai[wildcards.sample],
        germline_resource = gnomad_path.absolute()
    output:
        vcf = "mutect2/call/{sample}.vcf"
    message:
        "Calling variants on {wildcards.sample} with Mutect2"
    threads:
        2
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 20480)
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 115
        )
    params:
        extra = f"--germline-resource {str(gnomad_path.absolute())}"
    log:
        "logs/mutect2/call/{sample}.log"
    wrapper:
        f"{git}/bio/gatk/mutect"
