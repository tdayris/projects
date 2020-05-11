#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
This script splits a tsv line based on possible commas in some columns
"""

from typing import List

def iter_ann(ann: str) -> List[str]:
    """
    Return an iterator of fields within the annotation
    """
    if "," in ann:
        return ann.split(",")
    return [ann]


with open(snakemake.input[0], 'r') as infile, open(snakemake.output[0], 'w') as outfile:
    header = None
    for line in infile:
        if header is not None:
            content = {k: v for k, v in zip(header, line[:-1].split("\t"))}

            to_expand = {
                "effect": iter_ann(content["ANN[*].EFFECT"]),
                "gene": iter_ann(content["ANN[*].GENE"]),
                "gene_id": iter_ann(content["ANN[*].GENEID"]),
                "impact": iter_ann(content["ANN[*].IMPACT"]),
                "hgvs_p": iter_ann(content["ANN[*].HGVS_P"]),
                "hgvs_c": iter_ann(content["ANN[*].HGVS_C"]) ,
                "cdna_pos": iter_ann(content["ANN[*].CDNA_POS"]),
                "biotype": iter_ann(content["ANN[*].BIOTYPE"])
            }

            if not all(len(to_expand[k]) == len(to_expand["effect"]) for k in to_expand.keys()):
                raise ValueError("not same length")

            try:
                text = "\n".join([
                    "\t".join([
                        content["Sample_ID"],
                        content["Recurence"],
                        content["CHROM"],
                        content["POS"],
                        content["REF"],
                        content["ALT"],
                        content["FILTER"],
                        content["VAF"],
                        content["NAD"],
                        content["NRD"],
                        content["NDP"],
                        effect,
                        gene,
                        gene_id,
                        impact,
                        hgvs_p,
                        hgvs_c,
                        cdna_pos,
                        biotype,
                        content["dbNSFP_Interpro_domain"],
                        content["dbNSFP_1000Gp3_EUR_AC"],
                        content["dbNSFP_Uniprot_acc"],
                        content["dbNSFP_SIFT_pred"],
                        content["dbNSFP_MetaSVM_pred"],
                        content["dbNSFP_MutationTaster_pred"],
                        content["dbNSFP_MutationAssessor_pred"],
                        content["dbNSFP_PROVEAN_pred"],
                        content["dbNSFP_LRT_pred"],
                        content["dbNSFP_Polyphen2_HDIV_pred"],
                        content["dbNSFP_FATHMM_pred"],
                        content["GWASCAT_PUBMED_ID"],
                        content["MSigDb"],
                        content["dbNSFP_ExAC_NFE_AF"],
                        content["dbNSFP_ExAC_SAS_AF"],
                        content["dbNSFP_ExAC_Adj_AF"],
                        content["dbNSFP_ExAC_AFR_AF"],
                        content["dbNSFP_ExAC_AF"],
                        content["dbNSFP_ExAC_FIN_AF"],
                        content["dbNSFP_1000Gp3_EUR_AF"],
                        content["dbNSFP_ExAC_AMR_AF"],
                        content["dbNSFP_ExAC_EAS_AF"],
                        content["dbNSFP_ESP6500_EA_AF"],
                        content["dbNSFP_1000Gp3_AMR_AF"],
                        content["dbNSFP_1000Gp3_AF"],
                        content["dbNSFP_1000Gp3_EAS_AF"],
                        content["dbNSFP_1000Gp3_EUR_AF"],
                        content["dbNSFP_1000Gp3_AFR_AF"],
                        content["dbNSFP_1000Gp3_SAS_AF"],
                        content["MAX_AF"]
                    ])
                    for effect, gene, gene_id, impact,
                        hgvs_p, hgvs_c, cdna_pos, biotype
                    in zip(
                        to_expand["effect"],
                        to_expand["gene"],
                        to_expand["gene_id"],
                        to_expand["impact"],
                        to_expand["hgvs_p"],
                        to_expand["hgvs_c"],
                        to_expand["cdna_pos"],
                        to_expand["biotype"]
                    )
                ])
            except KeyError:
                print(line)
                raise
            else:
                outfile.write(f"{text}\n")

        else:
            header = line[:-1].split("\t")
            outfile.write(line)
