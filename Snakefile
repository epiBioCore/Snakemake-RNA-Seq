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
print(config)
def get_fastq(wildcards):
    return(units.loc[(wildcards.id,wildcards.unit),["fq1","fq2"]].dropna())


def get_replicate_fq(wildcards):
    input = expand([
            f"{config['outdir']}/Trimmed_fastq/{u.sample}_{u.unit}.{wildcards.read}.fq.gz"
            for u in units.loc[units["replicate"] == wildcards.replicate].itertuples()
        ])
    return(input)    




def get_fq(wildcards):
     if config["paired"]:
        return(expand(f"{config['outdir']}/Trimmed_fastq/{wildcards.sample}_{{read}}.fq.gz",read = [1,2]))
     else:
        return(f"{config['outdir']}/Trimmed_fastq/{wildcards.sample}_1.fq.gz")
    
wildcard_constraints:
    sample="|".join(replicates)


def get_strand(sample):
                report_file = [f for f in input if f.endswith(f"-{sample}.strandedness.txt")][0]
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



## fix rule all
rule all:
    input:
        #expand(config["outdir"] + "/Trimmed_fastq/{unit.sample}_{unit.unit}.1.fq.gz",unit=units.itertuples()),
        #expand(config["outdir"] + "/Trimmed_fastq/{unit.sample}_{unit.unit}.2.fq.gz",unit=units.itertuples()),
        #expand(config["outdir"] + "/FastQC/{unit.sample}_{unit.unit}_report.html",unit=units.itertuples()),
        #expand(config["outdir"] + "/Trimmed_fastq/{replicate}_{read}.fq.gz",replicate=replicates,read=[1,2])
        #expand(config["outdir"] + "/STAR/{sample}_Aligned.sortedByCoord.out.bam",sample=replicates)
        expand(config["outdir"] + "/STAR/{sample}.bam",sample=replicates),
        expand(config["outdir"] + "/Bigwigs/{sample}.bw",sample=replicates)

rule Fastqc:
    input: 
        get_fastq
                  
    output:
         html = config["outdir"] + "/FastQC/{id}_{unit}_report.html",
         zip = config["outdir"] + "/FastQC/{id}_{unit}_fastqc.zip"

    log:
        config["outdir"] + "/logs/fastqc/{id}_{unit}.log"

    wrapper:
        "0.65.0/bio/fastqc"



rule cutadapt_pe:
    input: 
        get_fastq

    output:
        fastq1 = config["outdir"] + "/Trimmed_fastq/{id}_{unit}.1.fq.gz",
        fastq2 = config["outdir"] + "/Trimmed_fastq/{id}_{unit}.2.fq.gz",
        qc = config["outdir"] + "/Trimmed_fastq/{id}_{unit}_qc.txt"

    log:
        config["outdir"] + "/logs/cutadapt/{id}_{unit}_cutadapt.out"
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
        read="1|2",
        replicate= "|".join(replicates)

    shell:
        "cat {input} > {output}"        



rule cutadapt_se:
    input: 
        get_fastq

    output:
        fastq1 = config["outdir"] + "/Trimmed_fastq/{sample}_{unit}.1.fq.gz",
        qc = config["outdir"] + "/Trimmed_fastq/{sample}_{unit}_qc.txt"

    log:
        config["outdir"] + "/logs/cutadapt/{sample}_{unit}_cutadapt.out"
    params:
        adapters = "-a " + config["cutadapt-se"]["r1_adapter"],
        others = config["cutadapt-se"]["params"]

    wrapper:
        "0.66.0/bio/cutadapt/se"


rule genome_generate:
    input:
        fasta = config["ref"],

    output:
        dir = directory("STARindex")    

    threads:
        10

    params:
        extra = "--sjdbGTFfile " + config["gtf"] + " --sjdbOverhang " + str(config["length"])

    wrapper:
        "0.66.0/bio/star/index"

rule STAR_mapping:
    input:
        fq= get_fq,
        genome = rules.genome_generate.output.dir

    output:
        bam = config["outdir"] + "/STAR/{sample}_Aligned.sortedByCoord.out.bam"

    params:
        outdir = config["outdir"] + "/STAR/{sample}_",   
    threads:
        config["star_mapping"]["threads"]

    log:
          config["outdir"] + "/logs/STAR/{sample}_STAR.out"
    shell:
          """
          STAR --runThreadN {threads} \
           --genomeDir {input.genome:} \
           --readFilesIn {input.fq} \
           --readFilesCommand gunzip -c \
           --runMode alignReads \
           --outSAMattributes All \
           --alignSJoverhangMin 8 \
           --alignSJDBoverhangMin 1 \
           --outFilterMismatchNmax 999 \
           --outFilterMismatchNoverLmax 0.04 \
           --alignIntronMin 20 \
           --alignIntronMax 1000000 \
           --alignMatesGapMax 1000000 \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --outFileNamePrefix {params.outdir} 2> {log}
        """  

rule index:
    input:
        bam = config["outdir"] + "/STAR/{sample}_Aligned.sortedByCoord.out.bam"

    output:
        bam = config["outdir"] + "/STAR/{sample}.bam"

    run:
        if config["paired"]:
            shell("samtools view -b -f 0x2 "+ "{input.bam}" + " > " + "{output.bam}")
            shell("samtools index " + "{output.bam}")
        else:
            shell("mv " + "{input.bam}" + " "+ "{output.bam}")   
            shell("samtools index " + "{output.bam}")

rule bigwigs:
    input:
        bam = config["outdir"] + "/STAR/{sample}.bam"

    output:
        bw = config["outdir"] + "/Bigwigs/{sample}.bw"

    log:
        config["outdir"] + "/logs/bamCoverage/{sample}.out"
    shell:"""
        bamCoverage -b {input.bam} --normalizeUsing CPM -of bigwig -o {output.bw}
        """

rule gtf_to_bed:
    input:
        config["gtf"]

    output:
        "resources/genome.bed" 

    shell:
        "echo"

rule infer_strand:
    input:
        bam = config["outdir"] + "/STAR/{sample}.bam",
        bed = "resources/genome.bed"

    output:
        txt = config["outdir"] + "/Strand/{sample}_strand.report.txt" 

    shell:
        "infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt}"       

rule cuffnorm:
    input:
        expand(config["outdir"] + "/STAR/{sample}.bam",sample=replicates)
    output:
        directory(config["outdir"] + "/Cuffnorm")

    params:
        gtf= config["gtf"], 
        lib= "first"       
    threads:
        config["cuffnorm"]["threads"]

    shell:""""
        cuffnorm -o $cuffnorm -p {threads} {params.gtf} {params.lib} {input}
        """
        
        