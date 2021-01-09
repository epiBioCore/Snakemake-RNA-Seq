from snakemake.utils import validate
import pandas as pd
import os.path
import glob

## inspired by https://github.com/snakemake-workflows/rna-seq-star-deseq2
configfile: "config.yaml"
units = pd.read_table(config["sample_sheet"], dtype=str)
if config["tech_reps"]:
    units.set_index(["Sample","tech_reps"], drop=False)
    #validate(units, schema="schemas/units.tech.reps.schema.yaml")

else:
    units.set_index("Sample")
    #validate(units, schema="schemas/units.schema.yaml")
    
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

#workdir: config["outdir"]

#samples = samples_df.index.get_level_values(0)
#units = samples_df.index.get_level_values(1)
#units = list(samples_df.unit.unique)

samples = units.Sample.tolist()
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





## fix rule all
rule all:
    input:
        #expand(config["outdir"] + "/Trimmed_fastq/{sample}.1.fq.gz",sample=samples),
        #expand(config["outdir"] + "/Trimmed_fastq/{sample}.2.fq.gz",sample=samples),
        expand(config["outdir"] + "/FastQC/{sample}_report.html",sample=samples)
        #expand(config["outdir"] + "/Trimmed_fastq/{sample}_{read}.fq.gz",sample=samples,read=[1,2])
        #expand(config["outdir"] + "/STAR/{sample}_Aligned.sortedByCoord.out.bam",sample=samples)
        #expand(config["outdir"] + "/STAR/{sample}.bam",sample=samples),
        #expand(config["outdir"] + "/Bigwigs/{sample}.bw",sample=samples),
        #directory(config["outdir"] + "/Cuffnorm"),



rule FastQC:
    input: 
        get_fastq
                  
    output:
         html = config["outdir"] + "/FastQC/{sample}_report.html",
         zip = config["outdir"] + "/FastQC/{sample}_fastqc.zip"

    log:
        config["outdir"] + "/logs/fastqc/{sample}.log"

    wrapper:
        "0.65.0/bio/fastqc"



rule trimmonatic_pe:
    input:
        get_fastq

    output:
        fastq1 = "{outdir}/Trimmed_fastq/{{sample}}_{fqext1}.{fqsuffix}".fomat(**config),
        fastq2 = "{outdir}/Trimmed_fastq/{{sample}}_{fqext2}.{fqsuffix}".fomat(**config),
        u1 = temp("{outdir}/Trimmed_fastq/{{sample}}_unpaired.1.fq.gz".format(**config)),
        u2 = temp("{outdir}/Trimmed_fastq/{{sample}}_unpaired.2.fq.gz".format(**config)),
        log = temp("{outdir}/Trimmed_fastq/{{sample}}_trimlog.txt.gz".format(**config))
        stat = "{outdir}/logs/Trimomatic/{{sample}}_trimmomatic.out".format(**config)

    benchmark:
        config["outdir"] + "/benchmarks/Trimomatic/{sample}_trimmomatic.out"

    params:
        partition = "talon-fat"

    resources:
        cpus = 16,
        time = "2:00:00",
        mem = "300G"

    shell: """
        trimmomatic PE -threads {resources.cpus} -phred33 -trimlog {output.log}  \
         {input} \
        {output.fastq1} {output.u1} \
         {output.fastq2} {output.u2} \
        ILLUMINACLIP:/home/danielle.perley/miniconda3/envs/atac-seq/share/trimmomatic-0.39-1/adapters/Truseq3-SE.fa:2:30:10:3:TRUE MINLEN:10 2> {output.stat}
    """

        
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

rule trimmonatic_se:
    input:
        get_fastq

    output:
        fastq1 = "{outdir}/Trimmed_fastq/{{sample}}_{fqext1}.{fqsuffix}".fomat(**config),
        log = temp("{outdir}/Trimmed_fastq/{{sample}}_trimlog.txt.gz".format(**config))
        stat = "{outdir}/logs/Trimomatic/{{sample}}_trimmomatic.out".format(**config)

  
    benchmark:
        "{outdir}/benchmarks/Trimomatic/{{sample}}_trimmomatic.benchmark".format(**config))

    params:
        partition = "talon-fat"

    resources:
        cpus = 16,
        time = "2:00:00",
        mem = "300G"

    shell: """
        trimmomatic SE -threads {resources.cpus} -phred33 -trimlog {output.log}  \
         {input} \
        {output.fastq1} \
        ILLUMINACLIP:/home/danielle.perley/miniconda3/envs/atac-seq/share/trimmomatic-0.39-1/adapters/Truseq3-SE.fa:2:30:10:3:TRUE MINLEN:10 2> {output.stat}
    """



if config["aligner"] == "STAR":

    rule genome_generate:
        input:
            fasta = config["ref"]

        output:
            dir = directory("STARindex")    

        resources:
            cpus = config["resources"]["star"]["cpus"],
            time = config["resources"]["star"]["time"],
            mem = config["resources"]["star"]["mem"]

        benchmark:
            config["outdir"] + "/benchmarks/STAR/genome_generate_bench.txt"

        log:
            config["outdir"] + "/logs/STAR/genome_generate_out.txt"
        params:
            gtf = config["gtf"],
            genome = config["ref"],
            len = config["length"],
            partition = "talon-fat"

        shell: """
                STAR --runThreadN {resources.cpus} \
                 --runMode genomeGenerate \
                --genomeDir {output.dir} \
                --genomeFastaFiles {params.genome} \
                --sjdbGTFfile {params.gtf} \
                --sjdbOverhang {params.len}
                """

    rule STAR_mapping:
        input:
            fq= get_trimmed,
            genome = rules.genome_generate.output.dir

        output:
            bam = config["outdir"] + "/Alignments/{sample}_sorted.bam",
            stat = config["outdir"] + "/Alignments/{sample}_Log.final.out"

        params:
            outdir = config["outdir"] + "/Alignments/{sample}_",   
        threads:
            config["star_mapping"]["threads"]

        benchmark:
            config["outdir"] + "/benchmarks/STAR/{sample}_STAR.bench.txt"

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

             mv {params.outdir}_Aligned.sortedByCoord.out.bam {output}
            """  

elif config["aligner"] == "hisat2":
    rule hisat:
        input:
            fq = get_fq,

        output:
            temp(bam = config["outdir"] + "/Alignments/{sample}.bam"),
            stat = config["outdir"] + "/Alignments/{sample}_align_stat.txt"
        params:
            input=(
                lambda wildcards, input: ["-U", input.fq]
                if config["Paired"] == False
                else ["-1", input.fq[0], "-2", input.fq[1]]),
            index = config["hisat2"]["index"],
            strand = config["hisat2"]["strand"],
            other = config["hisat2"]["other"],
            partition = "talon-fat"

        resources:
            cpus = config["resouces"]["hisat"]["cpus"],
            time = config["resouces"]["hisat"]["time"],
            mem = config["resouces"]["hisat"]["mem"],
   
        shell: """
         hisat2 -p {threads} --rna-strandness {params.strand} --new-summary {params.other} -x {params.index} {params.input} 2> {out.stat} | samtools view -hb - > {output.bam}
     
            """

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


rule cuffnorm:
    input:
        expand(config["outdir"] + "/Alignments/{sample}_sorted.bam",sample=samples)

    output:
        dir = directory(config["outdir"] + "/Cuffnorm"),
        table = directory(config["outdir"] + "/Cuffnorm/genes_fpkm.table")

    params:
        gtf= config["gtf"], 
        lib= config["cuffnorm"]["lib"],
        partition = "talon-fat"

    resources:
        cpus = config["resources"]["cuffnorm"]["cpus"],
        time = config["resources"]["cuffnorm"]["time"],
        mem = config["resources"]["cuffnorm"]["mem"],

    benchmark: 
        config["outdir"] + "/benchmarks/cuffnorm.bench.txt"
    
    log: 
        config["outdir"] + "/logs/cuffnorm.out"

    shell:"""
        cuffnorm -o {output.dir} -p {threads} --library-type {params.lib} {params.gtf} {input} 2> {log}
        """

rule featureCounts:
    input:
        expand(config["outdir"] + "/Alignment/{sample}_sorted.bam",sample=samples)

    output:
        counts = config["outdir"] + "/Counts/geneCounts.txt",
        stats = config["outdir"] + "/Counts/geneCounts.txt.summary"    
    
    log:
       config["outdir"] + "/logs/featureCounts.out" 
    benchmark:
       config["outdir"] + "/benchmarks/featureCounts.bench.txt" 
   

    params:
        gtf = config["gtf"],
        strand = config["featureCounts"]["strand"],
        other = config["featureCounts"]["other"]

    resources:
        cpus = config["resources"]["featureCounts"]["cpus"],
        time = config["resources"]["featureCounts"]["time"],
        mem = config["resources"]["featureCounts"]["mem"],


    shell: """
        featureCounts -a {params.gtf} -s {params.strand} {params.other} -o {output.counts} {input} 2> {log}
        """

rule prepfeatureCounts:
    input:
        counts = config["outdir"] + "/Counts/geneCounts.txt"

    output: 
        counts = config["outdir"] + "/Counts/geneCounts_for_DESEq2.csv",
        annot = config["outdir"] + "/Counts/gene_lengths.csv"

    script: "Script/prepfeatureCounts.R"    

      