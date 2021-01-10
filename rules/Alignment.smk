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
