def get_trimmed(wildcards):
    """
    Function that returns the reads for any aligner.
    """
    if config["paired"]:
        return sorted(expand("{outdir}/Trimmed_fastq/{{sample}}_{fqext}_Trimmed.{fqsuffix}", **config))
    return expand("{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}", **config)




if config["aligner"] == "STAR":

    rule genome_generate:
        input:
            fasta = config["ref"]

        output:
            dir = directory("{outdir}/STARindex".format(**config))    

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
                mkdir {output.dir}

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
            bam = "{outdir}/Alignments/{{sample}}_sorted.bam".format(**config),
            bai = "{outdir}/Alignments/{{sample}}_sorted.bam.bai".format(**config),
            stat = "{outdir}/Alignments/{{sample}}_Log.final.out".format(**config)

        params:
            outdir = "{outdir}/Alignments/{{sample}}_".format(**config), 
            partition = "talon-fat"  
        resources:
            cpus = config["resources"]["star"]["cpus"],
            time = config["resources"]["star"]["time"],
            mem = config["resources"]["star"]["mem"],

        benchmark:
            "{outdir}/benchmarks/STAR/{{sample}}_STAR.bench.txt".format(**config)

        log:
            "{outdir}/logs/STAR/{{sample}}_STAR.out".format(**config)
        shell:
            """
            STAR --runThreadN {resources.cpus} \
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

             mv {params.outdir}Aligned.sortedByCoord.out.bam {output.bam}
             samtools index {output.bam}
            """  

elif config["aligner"] == "hisat2":
    rule hisat:
        input:
            fq = get_trimmed,

        output:
            temp(bam = "{outdir}/Alignments/{{sample}}.bam".format(**config)),
            stat = "{outdir}/Alignments/{{sample}}_align_stat.txt".format(**config)
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
        bam = "{outdir}/Alignments/{{sample}}.bam".format(**config)

    output:
        bam = "{outdir}/Alignment/{{sample}}_sorted.bam".format(**config),
        bai = "{outdir}/Alignment/{{sample}}_sorted.bam.bai".format(**config)

    params:
        flag = "-f 0x2" if config["paired"] else "",
        partition = "talon-fat"

    benchmark:
        "{outdir}/benchmarks/samtools_sort_{{sample}}.bench.txt".format(**config)  

    resources:
        cpus = config["resources"]["samtools"]["cpus"],
        time = config["resources"]["samtools"]["time"],
        mem = config["resources"]["samtools"]["mem"],

    log:
        "{outdir}/logs/samtools_sort_{{sample}}.out".format(**config)   
    shell: """
        samtools view -hb {params.flag} {input.bam} | samtools sort -o {output.bam} - 
        samtools index {output.bam}
        """
