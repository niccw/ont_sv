#!/usr/bin/env python3
import argparse
import sys

def parse_out(out_path:str, ins_ls:list) -> list:
    """
    Parse the blast output (outfmt6) of SVIM ins/del sequence against te.fa
    """
    out_d = {}
    print(f"Parsing {out_path}...", file = sys.stderr)
    with open(out_path,"r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.strip().split("\t")
            sv_id = col[0].split("|")[0].split("_")[0] # each insertion seq is a uniqe record, all counted
            
            if sv_id not in ins_ls:
                continue
            
            # TE info
            te_superfamily = col[1].split("#")[1]
            align_len = int(col[3])
            if out_d.get(te_superfamily):
                out_d[te_superfamily] += align_len
            else:
                out_d[te_superfamily] = align_len
            
    return out_d

def parse_bed(bed_path:str,min_len:int = 10,min_support:int = 3):
    """
    Parse the bed file (from bedtools intersect sv_bed te_bed)
    Mark down (i.e. screen sv_len and sv_support) the INS to be screened in parse_out
    """
    bed_d = {}
    ins_ls = []
    print(f"Parsing {bed_path}...", file = sys.stderr)
    with open(bed_path,"r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.strip().split("\t")
            [sv_type,sv_len,sv_support,sv_id] = col[3].split(";")
            # skip all BND (breakpoint row); consistent with pairwise_compare_svim.py
            if sv_type == "BND" or sv_type == "INV":
                continue
            try:
                sv_len = abs(int(sv_len)) # length of deletion is negative
            except ValueError as ve:
                print("Invalid length value. Skip.\n{ve}", file = sys.stderr)

            sv_support = int(sv_support)

            if sv_len < min_len or sv_support < min_support:
                continue
            
            if sv_type == "INS":
                ins_ls.append(sv_id)

            # extract TE info and bp++
            te_superfamily = col[7].split(";")[1]
            if bed_d.get(te_superfamily):
                bed_d[te_superfamily] += int(col[10])
            else:
                bed_d[te_superfamily] = int(col[10])

    return [bed_d,set(ins_ls)]



def print_result(bed_d:dict, out_d:dict) -> None:
    print(f"# BED file: {args.bed}, OUT file: {args.out}", file = sys.stdout)
    # TE content calculated by intersect
    print(f"## TE content by intersect (exclude BND)", file = sys.stdout)
    for k,v in bed_d.items():
        print(f"{k}\t{v}", file = sys.stdout)

    # TE content calculated by blast
    print(f"## TE content by blast (for insertion sequence)", file = sys.stdout)
    for k,v in out_d.items():
        print(f"{k}\t{v}", file = sys.stdout)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type = str, dest= "out",help = "blast output (outfmt6) of SVIM ins/del sequence against te.fa")
    parser.add_argument("--bed", type = str, dest= "bed",help = "bed file (from bedtools intersect sv_bed te_bed)")
    parser.add_argument("--min_len", type = int, dest= "min_len",help = "default = 10", default = 10)
    parser.add_argument("--min_support", type = int, dest= "min_support",help = "default = 3", default = 3)
            
    args = parser.parse_args()

    bed_d, ins_ls = parse_bed(bed_path = args.bed, min_len = args.min_len, min_support = args.min_support)
    out_d = parse_out(out_path= args.out, ins_ls = ins_ls)
    print_result(bed_d = bed_d, out_d = out_d)
