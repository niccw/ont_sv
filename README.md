```
usage: pairwise_compare_svim_sv.py [-h] [-f FOLDER]
                                   [-k KEY_NAME [KEY_NAME ...]] [-w WINDOW]
                                   [-s SUPPORT] [-l SVLEN]

    vcf : *key_name*.vcf
    blast_out : *key_name*.out
    bedtools_inersect : *key_name*.bed
    out_name = f'te_annotated_{file_path_d["key"][i]}.bed'


optional arguments:
  -h, --help            show this help message and exit
  -f FOLDER             folder containning all files needed
  -k KEY_NAME [KEY_NAME ...]
                        sample key name
  -w WINDOW             window size add to INS (default = 20)
  -s SUPPORT            min. support (default = 3)
  -l SVLEN              min. SVLEN (default = 10)
  ```
