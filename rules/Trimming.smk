
rule trimmomatic_pe:
    input:
        r1 = "{fq_dir}/{{sample}}_{fqext1}.{fqsuffix}".format(**config),
        r2 = "{fq_dir}/{{sample}}_{fqext2}.{fqsuffix}".format(**config)


    output:
        fastq1 = "{outdir}/Trimmed_fastq/{{sample}}_{fqext1}_Trimmed.{fqsuffix}".format(**config),
        fastq2 = "{outdir}/Trimmed_fastq/{{sample}}_{fqext2}_Trimmed.{fqsuffix}".format(**config),
        u1 = temp("{outdir}/Trimmed_fastq/{{sample}}_unpaired.1.fq.gz".format(**config)),
        u2 = temp("{outdir}/Trimmed_fastq/{{sample}}_unpaired.2.fq.gz".format(**config)),
        log = temp("{outdir}/Trimmed_fastq/{{sample}}_trimlog.txt.gz".format(**config)),
        stat = "{outdir}/logs/Trimomatic/{{sample}}_trimmomatic.out".format(**config)

    benchmark:
        "{outdir}/benchmarks/Trimomatic/{{sample}}_trimmomatic.out".format(**config)

    params:
        adapter = config["trimmomatic-pe"]["adapters"],
        other = config["trimmomatic-pe"]["other"],
        partition = "talon-fat"

    resources:
        cpus = config["resources"]["trimmomatic-pe"]["cpus"],
        time = config["resources"]["trimmomatic-pe"]["time"],
        mem = config["resources"]["trimmomatic-pe"]["mem"]

    shell: """
        trimmomatic PE -threads {resources.cpus} -phred33 -trimlog {output.log}  \
         {input} \
        {output.fastq1} {output.u1} \
         {output.fastq2} {output.u2} \
        ILLUMINACLIP:{params.adapter}:2:30:10:3:TRUE {params.other} 2> {output.stat}
    """
rule trimmomatic_se:
    input:
        r1 = "{outdir}/{fq_dir}/{{sample}}_{fqext1}.{fqsuffix}".format(**config)

    output:
        fastq1 = "{outdir}/Trimmed_fastq/{{sample}}_Trimmed.{fqsuffix}".format(**config),
        log = temp("{outdir}/Trimmed_fastq/{{sample}}_trimlog.txt.gz".format(**config)),
        stat = "{outdir}/logs/Trimomatic/{{sample}}_trimmomatic.out".format(**config)

  
    benchmark:
        "{outdir}/benchmarks/Trimomatic/{{sample}}_trimmomatic.benchmark".format(**config)

    params:
        adapter = config["trimmomatic-se"]["adapters"],
        other = config["trimmomatic-se"]["other"],
        partition = "talon-fat"

    resources:
        cpus = config["resources"]["trimmomatic-se"]["cpus"],
        time = config["resources"]["trimmomatic-se"]["time"],
        mem = config["resources"]["trimmomatic-se"]["mem"]

    shell: """
        trimmomatic SE -threads {resources.cpus} -phred33 -trimlog {output.log}  \
         {input} \
        {output.fastq1} \
        ILLUMINACLIP:{params.adapter}:2:30:10:3:TRUE {params.other} 2> {output.stat}
    """
