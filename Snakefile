from snakemake.utils import validate
import pandas as pd
import os.path
import glob

## inspired by https://github.com/snakemake-workflows/rna-seq-star-deseq2
configfile: "config.yaml"
units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="schemas/units.schema.yaml")

#workdir: config["outdir"]

#samples = samples_df.index.get_level_values(0)
#units = samples_df.index.get_level_values(1)
#units = list(samples_df.unit.unique)
samples = list(units.index.get_level_values(0).unique())
names = units.descriptive_name.unique().tolist()
replicates = units.replicate.unique().tolist()

def get_fastq(wildcards):
    return(units.loc[(wildcards.sample,wildcards.unit),["fq1","fq2"]].dropna())

def is_single_end(sample,unit):
    return pd.isnull(units.loc[(sample,unit),"fq2"])

def get_replicate_fq(wildcards):
    input = expand([
            f"{config['outdir']}/Trimmed_fastq/{u.sample}_{u.unit}.{wildcards.read}.fq.gz"
            for u in units.loc[units["replicate"] == wildcards.replicate].itertuples()
        ])
    return(input)    




# def get_fq (wildcards):
#     if not is_single_end(**wildcards):
#         x
#     else:
#         y    

wildcard_constraints:
    sample = "|".join(samples)


## fix rule all
rule all:
    input:
        expand(config["outdir"] + "/Trimmed_fastq/{unit.sample}_{unit.unit}.1.fq.gz",unit=units.itertuples()),
        expand(config["outdir"] + "/Trimmed_fastq/{unit.sample}_{unit.unit}.2.fq.gz",unit=units.itertuples()),
        expand(config["outdir"] + "/FastQC/{unit.sample}_{unit.unit}_report.html",unit=units.itertuples()),
        expand(config["outdir"] + "/Trimmed_fastq/{replicate}_{read}.fq.gz",replicate=replicates,read=[1,2])



rule Fastqc:
    input: 
        get_fastq
                  
    output:
         html = config["outdir"] + "/FastQC/{sample}_{unit}_report.html",
         zip = config["outdir"] + "/FastQC/{sample}_{unit}_fastqc.zip"

    log:
        config["outdir"] + "/logs/fastqc/{sample}_{unit}.log"

    wrapper:
        "0.65.0/bio/fastqc"



rule cutadapt_pe:
    input: 
        get_fastq

    output:
        fastq1 = config["outdir"] + "/Trimmed_fastq/{sample}_{unit}.1.fq.gz",
        fastq2 = config["outdir"] + "/Trimmed_fastq/{sample}_{unit}.2.fq.gz",
        qc = config["outdir"] + "/Trimmed_fastq/{sample}_{unit}_qc.txt"

    log:
        config["outdir"] + "/logs/cutadapt/{sample}_{unit}_cutadapt.out"
    params:
        adapters = "-a " + config["cutadapt-pe"]["r1_adapter"] + " -g " + config["cutadapt-pe"]["r2_adapter"],
        others = config["cutadapt-pe"]["params"]

    wrapper:
        "0.66.0/bio/cutadapt/pe"
        
rule merge_tech_reps:
    input:
        get_replicate_fq

    output:
        config["outdir"] +"/Trimmed_fastq/{replicate}_{read}.fq.gz"

    wildcard_constraints:
        read="1|2"

    shell:
        "cat {input} > {output}"        



# rule trim_galore-se:
#     input: 
#         get_fastq

#     output:
#         r1 = config["outdir"] + "/Trimmed_fastq/{sample}_{unit}_val_1.fq.gz",

#     log:
#         config["outdir"] + "/logs/trim/{sample}_{unit}_trim_galore.out"
#     params:
#         o = config["outdir"] + "/Trimmed_fastq",
#         extra = config["trim_galore-se"]["params"]

#     shell:
#         "trim_galore {params.extra} -o {params.o} -j 8 {input}"


 
# rule:
#     input:
#         r1 = "Grove_results/trim/Trim_galore/{sample}_R1_val_1.fq.gz",
#         r2 = "Grove_results/trim/Trim_galore/{sample}_R2_val_2.fq.gz"

#     output:
#         bam = "Grove_results/Alignment/{sample}_trim_galore_Aligned.sortedByCoord.out.bam"

#     params:
#         outdir = "Grove_results/Alignment/{sample}_trim_galore_"   

#     log:
#         "Grove_results/logs/Alignment/{sample}_trim_galore_STAR.out"
#     shell:
#          """
#          STAR --runThreadN 10 \
#           --genomeDir ../RNA-Seq/STAR1index \
#           --readFilesIn {input.r1} {input.r2} \
#           --readFilesCommand gunzip -c \
#           --runMode alignReads \
#           --outSAMattributes All \
#           --alignSJoverhangMin 8 \
#           --alignSJDBoverhangMin 1 \
#           --outFilterMismatchNmax 999 \
#           --outFilterMismatchNoverLmax 0.04 \
#           --alignIntronMin 20 \
#           --alignIntronMax 1000000 \
#           --alignMatesGapMax 1000000 \
#           --outSAMstrandField intronMotif \
#           --outSAMtype BAM SortedByCoordinate \
#           --outSAMunmapped Within KeepPairs \
#           --outFileNamePrefix {params.outdir}
#        """  
#