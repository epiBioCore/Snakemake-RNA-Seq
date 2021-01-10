from snakemake.utils import validate
import pandas as pd
import os.path
import glob

## inspired by https://github.com/snakemake-workflows/rna-seq-star-deseq2
configfile: "config.yaml"
units = pd.read_table(config["sample_sheet"], dtype=str).set_index("Sample")
    #validate(units, schema="schemas/units.schema.yaml")
    
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

#workdir: config["outdir"]

#samples = samples_df.index.get_level_values(0)
#units = samples_df.index.get_level_values(1)
#units = list(samples_df.unit.unique)

samples = units.index
print(samples)
## coreNumbers must be unique
assert len(samples) == len(set(samples)),"Duplicate sample names detected. Check your sample sheet."

#names = units.descriptive_name.unique().tolist()
#replicates = units.replicate.unique().tolist()
print(config)
def get_fastq(wildcards):
    if config["paired"]:
        r1 = "{fq_dir}/{{sample}}_{fqext1}.{fqsuffix}".format(**config)
        r2 = "{fq_dir}/{{sample}}_{fqext2}.{fqsuffix}".format(**config)
        return([r1,r2])
    else:
        r1 = "{fq_dir}/{{sample}}_{fqext1}.{fqsuffix}".format(**config),
        return(r1)


def get_replicate_fq(wildcards):
    input = expand([
            f"{config['outdir']}/Trimmed_fastq/{u.sample}_{u.unit}.{wildcards.read}.fq.gz"
            for u in units.loc[units["replicate"] == wildcards.replicate].itertuples()
        ])
    return(input)    




def get_trimmed(wildcards):
     if config["paired"]:
        return(expand(f"{config['outdir']}/Trimmed_fastq/{wildcards.sample}_{{read}}.fq.gz",read = [1,2]))
     else:
        return(f"{config['outdir']}/Trimmed_fastq/{wildcards.sample}_1.fq.gz")
    
wildcard_constraints:
    sample="|".join(samples)

include: "rules/Alignment.smk"
include: "rules/Trimming.smk"
include: "rules/QC.smk"
include: "rules/Quantification.smk"
include: "rules/merge_replicates.smk"




## fix rule all
rule all:
    input:
        expand(config["outdir"] + "/Trimmed_fastq/{sample}_R1_Trimmed.fastq.gz",sample=samples),
        expand(config["outdir"] + "/Trimmed_fastq/{sample}_R2_Trimmed.fastq.gz",sample=samples),
        expand(config["outdir"] + "/FastQC/{sample}_report.html",sample=samples)
        #expand(config["outdir"] + "/Trimmed_fastq/{sample}_{read}.fq.gz",sample=samples,read=[1,2])
        #expand(config["outdir"] + "/STAR/{sample}_Aligned.sortedByCoord.out.bam",sample=samples)
        #expand(config["outdir"] + "/STAR/{sample}.bam",sample=samples),
        #expand(config["outdir"] + "/Bigwigs/{sample}.bw",sample=samples),
        #directory(config["outdir"] + "/Cuffnorm"),








        
rule merge_tech_reps:
    input:
        get_replicate_fq

    output:
        config["outdir"] +"/Trimmed_fastq/{replicate}_{read}.fq.gz"

    wildcard_constraints:
        read="1|2",
        replicate= "|".join(samples)

    shell:
        "cat {input} > {output}"        




rule sort:
    input:
        bam = config["outdir"] + "/Alignments/{sample}.bam"

    output:
        bam = config["outdir"] + "/Alignment/{sample}_sorted.bam"

    params:
        flag = "-f 0x2" if config["paired"] else "",
        partition = "talon-fat"

    benchmark:
        config["outdir"] + "/benchmarks/samtools_sort_{sample}.bench.txt"  

    resources:
        cpus = config["resources"]["samtools"]["cpus"],
        time = config["resources"]["samtools"]["time"],
        mem = config["resources"]["samtools"]["mem"],

    log:
        config["outdir"] + "/logs/samtools_sort_{sample}.out"   
    shell: """
        samtools view -hb {params.flag} {input.bam} | samtools sort -o {output.bam} - 
        samtools index {output.bam}
        """


rule bigwigs:
    input:
        bam = config["outdir"] + "/Alignments/{sample}_sorted.bam"

    output:
        bw = config["outdir"] + "/Bigwigs/{sample}.bw"

    log:
        config["outdir"] + "/logs/bamCoverage/{sample}.out"
    
    benchmark:
        config["outdir"] + "/benchmarks/bamCoverage/{sample}.bench.txt"
    resources:
        cpus = config["resources"]["bamCov"]["cpus"],
        time = config["resources"]["bamCov"]["time"],
        mem = config["resources"]["bamCov"]["mem"],


    shell:"""
        bamCoverage -b {input.bam} --normalizeUsing CPM -of bigwig -o {output.bw}
        """

