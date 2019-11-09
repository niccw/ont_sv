#!/usr/bin/env python3
import sys
import re

"""
Usage: ngmlr_svim_vcf2bed.py [svim.vcf]
Output:
    BED format to stdout
    The start position of all SVs are col[1] (VCF and BED are both 0 base)
    For BND and INS, the end position is 1 + start position (as they have no
    length); else are END tag
    on col[7] info field.
"""

def ngmlr_vcf2bed(path:str):
    with open(path,"r") as f:
        for line in f:
            if line.startswith("##") or line.startswith("#"):
                continue
            col = line.strip().split("\t")
            
            sv_type = re.search(r"(?<=SVTYPE=)\w+", col[7])[0]
            sv_len_m = re.search(r"(?<=SVLEN=)(.*?)(?=;)", col[7])
           
           # Parse info field (col[7])
            try:
                infos = parse_info(col[7])
            except Exception as e:
                print(col)
                sys.exit("cant parse col7 info")

            if sv_len_m is None:
                sv_len = "NA"
            else:
                sv_len = sv_len_m[0]
            
            if infos["SVTYPE"] in ["BND","INS"]:
                # BND and INS do not have corresponding coord. on ref. genome
                sv_end = int(col[1]) + 1
            else:
                sv_end = int(infos["END"])

            # output bed 4 col formatH            # chrom start   end
            # [sv_type;sv_len;sv_qual;sv-id]
            print(f'{col[0]}\t{col[1]}\t{sv_end}\t{sv_type};{sv_len};{infos["SUPPORT"]};{col[2]}')

def parse_info(s:str)->dict:
    d = {}
    d["SVTYPE"] = re.search(r'(?<=SVTYPE=)\w+',s).group(0)
    d["SUPPORT"] = re.search(r'(?<=SUPPORT=)\d+',s).group(0)
    if d["SVTYPE"] in ["BND"]:
        return d
    d["END"] = re.search(r'(?<=END=)\d+',s).group(0)
    if d["SVTYPE"] in ["INV"]:
        return d
    d["SVLEN"] = re.search(r'(?<=SVLEN=)(.*?)(?=;)',s).group(0)
    d["READS"] = re.search(r'(?<=READS=)(.*?)(?=$)',s).group(0).split(",")
    if d["SVTYPE"] == "INS":
        d["SEQS"] = re.search(r'(?<=SEQS=)(.*?)(?=;)',s).group(0).split(",")
    return d

if __name__ == "__main__":
    ngmlr_vcf2bed(sys.argv[1])
