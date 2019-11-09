from pybedtools import BedTool
import sys

# test pybedtools

f1 = BedTool("te_annotated_105.bed")

# This works :)
result = f1.intersect(b=["te_annotated_3A.bed", "te_annotated_6A.bed"],wao=True,names=["SampleA","SampleB"])

print("hi")