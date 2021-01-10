ruleorder: trimmonatic_pe > trimmonatic_se

rule trimmonatic_pe:
    input:
        get_fastq

    output:
        fastq1 = "{outdir}/Trimmed_fastq/{{sample}}_{fqext1}_Trimmed.{fqsuffix}".format(**config),
        fastq2 = "{outdir}/Trimmed_fastq/{{sample}}_{fqext2}_Trimmed.{fqsuffix}".format(**config),
        u1 = temp("{outdir}/Trimmed_fastq/{{sample}}_unpaired.1.fq.gz".format(**config)),
        u2 = temp("{outdir}/Trimmed_fastq/{{sample}}_unpaired.2.fq.gz".format(**config)),
        log = temp("{outdir}/Trimmed_fastq/{{sample}}_trimlog.txt.gz".format(**config)),
        stat = "{outdir}/logs/Trimomatic/{{sample}}_trimmomatic.out".format(**config)

    benchmark:
        config["outdir"] + "/benchmarks/Trimomatic/{sample}_trimmomatic.out"

    params:
        adapter = config["trimmomatic-pe"]["adapters"],
        other = config["trimmonatic-pe"]["other"],
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
        ILLUMINACLIP:{params.adapter}:2:30:10:3:TRUE {params.other} 2> {output.stat}
    """
rule trimmonatic_se:
    input:
        get_fastq

    output:
        fastq1 = "{outdir}/Trimmed_fastq/{{sample}}_Trimmed.{fqsuffix}".format(**config),
        log = temp("{outdir}/Trimmed_fastq/{{sample}}_trimlog.txt.gz".format(**config)),
        stat = "{outdir}/logs/Trimomatic/{{sample}}_trimmomatic.out".format(**config)

  
    benchmark:
        "{outdir}/benchmarks/Trimomatic/{{sample}}_trimmomatic.benchmark".format(**config)

    params:
        adapter = config["trimmomatic-se"]["adapters"],
        other = config["trimmonatic-se"]["other"],
        partition = "talon-fat"

    resources:
        cpus = 16,
        time = "2:00:00",
        mem = "300G"

    shell: """
        trimmomatic SE -threads {resources.cpus} -phred33 -trimlog {output.log}  \
         {input} \
        {output.fastq1} \
        ILLUMINACLIP:{params.adapter}:2:30:10:3:TRUE {params.other} 2> {output.stat}
    """
