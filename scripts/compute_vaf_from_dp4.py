#!/usr/bin/python3

import re

info_written = False
regex = r"DP4=[^;]+;"
new_header = ('##INFO=<ID=AF,Number=1,Type=Float,'
              'Description="Ratio between reads supporting alternative, '
              'and reads supporting reference">')


with open(snakemake.input["vcf"], r) as infile, open(snakemake.output["vcf"] as outfile):
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            if not line.startswith("##") and info_written is False:
                outfile.write(new_header)
                info_written = True
            continue

        chomp = line.split("\t")
        nf, fr, af, ar = re.match(regex, chomp[8]).group().split(",")
        chomp[8] += f";AF={(nf+fr)/(af+ar)}"
        outfile.write("\t".join(chomp))
