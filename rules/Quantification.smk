
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

      