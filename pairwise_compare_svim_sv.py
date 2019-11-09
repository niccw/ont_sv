import sys
import os
import re
import argparse
import itertools
from pybedtools import BedTool

"""
This script ignore breakpoint SV (BND)
"""

# check files
def check_file_path(folder:str,key_names:list) -> dict:
    d = {
        "key" : [],
        "vcf" : [],
        "blast" :[],
        "intersect" : []}

    file_list = [r + "/" + i for r,_,x in os.walk(folder) for i in x]
    #file_list = [y for x in file_list for y in x] # flatten the list
    #_,_,file_list = os.walk(folder)
    for key in key_names:
        d["key"].append(key)
        # re.compile(rf".*{key}*vcf")
        # list(filter(lambda x: re.search(rf"{key}.*vcf", x), file_list))
        vcf_re = list(filter(lambda x: re.search(rf"{key}.*vcf", x), file_list))
        out_re = list(filter(lambda x: re.search(rf"{key}.*out", x), file_list))
        bed_re = list(filter(lambda x: re.search(rf"{key}.*bed", x), file_list))
        if vcf_re:d["vcf"].append(vcf_re[0])
        if out_re:
            d["blast"].append(out_re[0])
        else:
            d["blast"].append(None)
        if bed_re:
            d["intersect"].append(bed_re[0])
        else:
            d["intersect"].append(None)

    if len(d["vcf"]) != len(d["key"]): # must found all VCF
        sys.exit(f'Cannot find all VCF.\n Keys entered: {",".join(key_names)}\nVCF found: {",".d["vcf"]}')
    return d

# read TE info from out & bed
def vcf2bed_annotateTE(file_path_d:dict, window:int = 20) -> None:
    """
    window is added for better filtering of shared INS across sample. (the insert position varies a little bit)
    """
    blast_d = None
    intersect_d = None

    for i in range(len(file_path_d["key"])):
        # for each experiment/condition, read blast and bedtools output
        if file_path_d["blast"][i] is not None:
            blast_d = te_info2map(file_path_d["blast"][i],"blast")
        else:
            print(f'file_path_d["key"][i]: no blast.out available, skip.')
        if file_path_d["intersect"][i] is not None:
            intersect_d = te_info2map(file_path_d["intersect"][i],"intersect")
        else:
            print(f'file_path_d["key"][i]: no bedtools intersect.bed available, skip.')

        out_name = f'te_annotated_{file_path_d["key"][i]}_INSpad{window}.bed'
        if os.path.exists(out_name):
            q = input(f"te_annotated_{out_name}.vcf already exist, rewrite it? (Y/N)")
            if q.capitalize() == "N":
                sys.exit(0)

        print(f'Open fh on {out_name}. Convert VCF to BED (read comments in script for details), subset of INFO parse to NAME (col 4) field.',file=sys.stderr)

        with open(file_path_d["vcf"][i],"r") as f, open(f'{out_dir}/{out_name}', "w") as o:
            line_count = 1
            for line in f:
                line_count += 1
                if line.startswith("##") or line.startswith("#"):
                    continue
                col = line.strip().split("\t")
                try:
                    infos = parse_info(col[7])
                except Exception as e:
                    print(f"{line_count}: Cannot parse info field.\n{line}\n{e}")
                
                sv_chr = col[0]
                sv_start = int(col[1])          
                sv_end = int(sv_start) + 1 if "END" not in infos else int(infos["END"]) # if missing END (i.e. BND) use start + 1
                sv_id = col[2]

                name = f'ID={sv_id};SVTYPE={infos["SVTYPE"]};SUPPORT={infos["SUPPORT"]}'
                if "SVLEN" in infos:
                    name += f';SVLEN={infos["SVLEN"]}'

                # chr start end name{ID;SVTYPE;SUPPORT;SVLEN;BLAST_TE (sep=,);INTERSECT_TE(sep=,)}
                if infos["SVTYPE"] == "INS":
                    sv_start = sv_start - 10 if (sv_start - 10) > 0 else 0
                    sv_end = sv_end + 10  # there is chance that sv_end larger than chr length, but should be rare and we can filter this later
                    if blast_d is not None:
                        if sv_id in blast_d:
                            name += f';BLAST_TE={blast_d[sv_id]}'
                if intersect_d is not None:
                    if sv_id in intersect_d:
                        name += f';INTERSECT_TE={intersect_d[sv_id]}'
                  
                # write to out_file
                # if missing END (i.e. BND) use start + 1
                o.write(f'{sv_chr}\t{sv_start}\t{sv_end}\t{name}\n')
        print(f'Finish writing {out_name}. Close fh.',file=sys.stderr)
                            
# bed/blast annotation to map
# return dict[sv_id:str] = te:str
def te_info2map(file_path:str,file_type:str) -> dict:
    d = {}
    with open(file_path,"r") as f:
        for line in f:
            col = line.strip().split("\t")
            # sv_id, sv_chr
            if file_type == "blast":
                sv_id = col[0].split("|")[0]
                te = col[1]
                # take care about INS, add ALL TE match to one key-value, as string (sep=,)
                if "_" in sv_id: # e.g. svim.INS.57_2
                    sv_id = sv_id.split("_")[0]
            elif file_type == "intersect":
                te = col[7]
                sv_id = col[3].split(";")[3]

            # fix te seperator ";" -> "#"
            if ";" in te:
                te = re.sub(";","#",te)
            
            if sv_id not in d:
                d[sv_id] = te
            else:
                if te in d[sv_id]:
                    continue
                else:
                    d[sv_id] = d[sv_id] + "," + te
    return d

# pairwise intersect of bed
def bed_comp(file_path_d:dict):
    # check and print all annotated BED which are going to be processed.
    file_path_d["annotated_bed"] = []
    file_list = [i for _,_,x in os.walk(f'{out_dir}') for i in x]
    for i in range(len(file_path_d["key"])):
        m = list(filter(lambda x: re.search(rf'te_annotated_{file_path_d["key"][i]}.*.bed', x), file_list))[0]
        if not m:
            sys.exit(f'Cannot locate te-annotated bed for {i}')
        else:
            file_path_d["annotated_bed"].append(f'{out_dir}/{m}')
        """
        annotated_bed = f'te_annotated_{file_path_d["key"][i]}_INSpad{window}.bed'
        file_path_d["annotated_bed"].append(annotated_bed)
        if not os.path.exists(f'{out_dir}/{annotated_bed}'):
            sys.exit(f'Cannot locate {annotated_bed}.')
        """
    
    bed_out_dir = "bed_intersect_output"
    if not os.path.exists(f'{out_dir}/{bed_out_dir}'):
        os.mkdir(f'{out_dir}/{bed_out_dir}')

    for a_bed, key in zip(file_path_d["annotated_bed"], file_path_d["key"]):
        beds_to_compare = [i for i in file_path_d["annotated_bed"] if i != a_bed]
        keys_to_compare = [i for i in file_path_d["key"] if i != key]
        bed1 = BedTool(a_bed)

        #result = f1.intersect(b=["te_annotated_3A.bed", "te_annotated_6A.bed"],wao=True,names=["SampleA","SampleB"])
        bed_wao_multi_intersect  = bed1.intersect(b=beds_to_compare, names=keys_to_compare,wao=True)
        bed_wao_multi_intersect.moveto(f'{out_dir}/{bed_out_dir}/mutiintersect_{key}.bed')


# generate report
def check_overlap(key_list:list, support:int = 3, sv_len:int = 10):
    sorted_key_list = sorted(key_list)
    # f'{out_dir}/mutiintersect_{key}.bed
    # check if out dir exist
    intersect_bed_dir = f'{out_dir}/bed_intersect_output'
    if not os.path.exists(intersect_bed_dir):
        sys.exit("bed_intersect_output/ not exist, previous steps probably went wrong.")
    
    intersect_bed_list = [x for _,_,i in os.walk(intersect_bed_dir) for x in i]

    for bed in intersect_bed_list:
        subject_key = re.search(r'(?<=mutiintersect_).+(?=\.bed)',bed)[0]
        if subject_key not in key_list:
            print(f'Found {bed}, but thsi sample is not in input sampel key name. Skip.')

    """
    # prepare a dict to marked the shared and unique SV NUMBER
     key_combination_ls = []
    for i in range(len(key_list)):
        key_combination_ls.extend(list(itertools.combinations(key_list, i+1)))
    key_combination_ls = [tuple(sorted(t)) for t in key_combination_ls]
    key_combination_cntd = dict(zip(key_combination_ls,[0] * len(key_combination_ls)))
    """

    # Mark down all combination for EACH SV
    cnt_d = {} # cnt sv intersect cnt_d[<key>][<sv_id>] = [...](0/1 for each pair of possible intersect)
    
    for bed in intersect_bed_list:
        subject_key = re.search(r'(?<=mutiintersect_).+(?=\.bed)',bed)[0]
        keys_to_compare = [i for i in key_list if i != subject_key]
        with open(f'{intersect_bed_dir}/{bed}',"r") as f, open(f'{out_dir}/unique_sv_{subject_key}.id',"w") as o:
            for line in f:
                col = line.strip().split("\t")
                sv_id = re.search(r'(?<=ID=).*?(?=;)',col[3])[0]
                
                subject_support = re.search(r'(?<=SUPPORT=)\d+',col[3]).group(0)
                subject_sv_type = re.search(r'(?<=SVTYPE=)\w+',col[3]).group(0)
                try:
                 subject_len = re.search(r'(?<=SVLEN=)(.*?)(?=;)|(?<=SVLEN=)(.*)(?=$)',col[3]).group(0)
                except AttributeError:
                    subject_len = 0

                if int(subject_support) >= support and int(subject_len) >= sv_len:
                    # build subject_sample|sv_id|SUPPORT
                    subject_key_sv_id = f'{subject_key}|{sv_id}|{subject_support}|{subject_len}'
                    if subject_key_sv_id not in cnt_d:
                        cnt_d[subject_key_sv_id] = []
                    
                    # col 5 is the key to compare, if "." -> not match to all other samples
                    if col[4] == ".":
                        o.write(f'{sv_id}\n')
                    else:
                        # maybe other_d[sorted tuple(sv_id)] = min.support among all
                        
                        # check target SUPPORT
                        target_support = re.search(r'(?<=SUPPORT=)\d+',col[8]).group(0)
                        target_sv_type = re.search(r'(?<=SVTYPE=)\w+',col[8]).group(0)
                        try:
                            target_len = re.search(r'(?<=SVLEN=)(.*?)(?=;)|(?<=SVLEN=)(.*)(?=$)',col[8]).group(0)
                        except AttributeError:
                            target_len = 0

                        if int(target_support) >= support and subject_sv_type == target_sv_type and int(target_len) >= sv_len:
                            target_sv_id = re.search(r'(?<=ID=).*?(?=;)',col[8])[0]
                            target_key_sv_id = f'{col[4]}|{target_sv_id}|{target_support}|{target_len}'
                            cnt_d[subject_key_sv_id].append(target_key_sv_id)

    # flatten the dict
    # TODO: add TE info
    #out_name = f'te_annotated_{file_path_d["key"][i]}.bed'
    te_annotation_d = {}
    for bed in [f for f in os.listdir(f'{out_dir}') if re.match(r'te_annotated.*',f)]:
        sample_key = re.search(r'(?<=te_annotated_).*?(?=_)',bed)[0]
        if sample_key not in te_annotation_d:
            te_annotation_d[sample_key] = {}  # te_annotation['105'][svim.DEL.1] = [<blast_te_anno>, <intersect_te_anno>]
        with open(f'{out_dir}/{bed}',"r") as f:
            for line in f:
                col = line.strip().split("\t")
                svim_id = re.search(r'(?<=ID=).*?(?=;)',col[3])[0]
                blast_te = ""
                intersect_te = ""
                try:
                    blast_te = re.search(r'(?<=BLAST_TE=).*?(?=;)|(?<=BLAST_TE=).*?(?=$)',col[3]).group(0)
                except AttributeError:
                    pass
                try:
                    intersect_te = re.search(r'(?<=INTERSECT_TE=).*?(?=;)|(?<=INTERSECT_TE=).*?(?=$)',col[3]).group(0)
                except AttributeError:
                    pass
                te_annotation_d[sample_key][svim_id] = [blast_te,intersect_te]


    sv_set = {}
    for k,v in cnt_d.items(): # v is the count list for each sv
        # build value for sv_set
        sample_key_ls = []
        sample_key_ls.append(k.split("|")[0])
        for i in v:
            sample_key_ls.append(i.split("|")[0])
        sample_key_ls = sorted(set(sample_key_ls))

        # build key(tuple) for sv_set
        k = [k]
        k.extend(v)
        k.sort()
        k = tuple(k)

        if k not in sv_set:
            # add TE info
            blast_te = []
            intersect_te = []
            for sv_id in k:
                [sample_key,svim_id] = sv_id.split("|")[0:2]
                blast_te.append(te_annotation_d[sample_key][svim_id][0])
                intersect_te.append(te_annotation_d[sample_key][svim_id][1])
                
            sv_set[k] = [sample_key_ls,blast_te,intersect_te]

    # report Summary of SV using key_combination_cntd
    with open(f'{out_dir}/key_combination_sv_cnt.tsv',"w") as o:
        o.write(f'##Directory: {args.folder}; Keys: {args.key_name}; INS padding size: {args.window}; min. support: {args.support}; min. svlen: {args.svlen}\n')
        for k,v in sv_set.items():
            o.write(f'{",".join(v[0])}\t{",".join(k)}\t{";".join(v[1])}\t{";".join(v[2])}\n')
        

#parse the key_combination_sv_cnt.tsv and report the number of SV for each combination
def final_summary(cnt_output:str = "key_combination_sv_cnt.tsv"):
    key_combination_d = {} # use tupple as key i.e. ('105', '3A')
    with open(f'{out_dir}/{cnt_output}',"r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            [comb,sv_id] = line.strip().split("\t")[0:2]
            comb = tuple(comb.split(","))
            sv_id = sv_id.split(",")
            sv_type = ""
            for i in sv_id:
                [sample,svim_id,support,sv_len] = i.split("|")
                i_sv_type = svim_id.split(".")[1]
                if i_sv_type != sv_type: sv_type = i_sv_type
            if sv_type not in key_combination_d:
                key_combination_d[sv_type] = {}
            if comb not in key_combination_d[sv_type]:
                key_combination_d[sv_type][comb] = 1
            else:
                key_combination_d[sv_type][comb] += 1

        with open(f'{out_dir}/summary.txt',"w") as f:
            for sv_type, comb_dict in key_combination_d.items():
                f.write(f'##{sv_type}\n')
                for k,v in comb_dict.items():
                    f.write(f'{k}\t{v}\n')

        
## helper function ##
def parse_info(s:str) -> dict:
    """
    # parse vcf col[7] info field
    BND: SVTYPE, SUPPORT
    INV: SVTYPE, SUPPORT, END
    OTHER (except INS): SVTYPE, SUPPORT, END, SVLEN, READS
    INS: SVTYPE, SUPPORT, END, SVLEN, READS, SEQS
    """
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
    # argsparse
    help_text = """
    vcf : *key_name*.vcf
    blast_out : *key_name*.out
    bedtools_inersect : *key_name*.bed
    out_name = f'te_annotated_{file_path_d["key"][i]}.bed'
    """
    parser = argparse.ArgumentParser(description=help_text,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", type = str, dest= "folder",help = "folder containning all files needed", default= "")
    parser.add_argument("-k", dest= "key_name",help = "sample key name", nargs= "+")
    parser.add_argument("-w", type = int, dest= "window",help = "window size add to INS (default = 20)", default = 20)
    parser.add_argument("-s", type = int, dest= "support",help = "min. support (default = 3)", default = 3)
    parser.add_argument("-l", type = int, dest= "svlen",help = "min. SVLEN (default = 10)", default = 10)
    args = parser.parse_args()

    out_dir = "compare_svim_output"
    # automatically check if all necessary files exist
    file_path_d = check_file_path(folder=args.folder,key_names=args.key_name)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    # parse the VCF to BED, adding TE annotation from blast and bedtool intersect out, also add padding for INS
    vcf2bed_annotateTE(file_path_d=file_path_d, window = args.window)
    # intersect the BED using pybedtool
    bed_comp(file_path_d = file_path_d)
    # parse and summarise BED intersect output in <group>\t<[svid]>
    check_overlap(key_list = args.key_name, support = args.support, sv_len = args.svlen)
    # Do the counting
    final_summary()

        