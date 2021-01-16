import pandas as pd
import os.path
import re
import csv

def get_strand(report_file):
    with open(report_file) as report:
        fail_val = fwd_val = 0
        for line in report:
            if line.startswith("Fraction of reads failed"):
                fail_val = float(line.strip().split(": ")[1])
            elif line.startswith(("""Fraction of reads explained by "1++""",
                             """Fraction of reads explained by "++""")):
                fwd_val = float(line.strip().split(": ")[1])
    if fwd_val > 0.6:
        return "forward"
    elif 1 - (fwd_val + fail_val) > 0.6:
        return "reverse"
    else:
        return "no"

strand = {}

#files = ["/lower_bay/home/danielle.perley/Snakemake-RNA-Seq/RNA-Seq_results/strandness/control_1_strand.report.txt","/lower_bay/home/danielle.perley/Snakemake-RNA-Seq/RNA-Seq_results/strandness/Treat_1_strand.report.txt"]
for file in snakemake.input:
    base = os.path.basename(file)
    sample = re.sub("_strand.report.txt","",base)

    strand[sample] = get_strand(file)

strandedness = pd.DataFrame.from_dict(strand, orient = 'index',columns = ["Strandedness"])


strandedness.to_csv(snakemake.output[0],sep="\t")
 
