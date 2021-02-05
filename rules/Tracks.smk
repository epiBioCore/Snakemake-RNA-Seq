rule bigwigs:
    input:
        bam = "{outdir}/Alignments/{{sample}}_sorted.bam".format(**config)

    output:
        bw = "{outdir}/Bigwigs/{{sample}}.bw".format(**config)

    log:
        "{outdir}/logs/bamCoverage/{{sample}}.out".format(**config)
    
    benchmark:
        "{outdir}/benchmarks/bamCoverage/{{sample}}.bench.txt".format(**config)
    
    params:
        partition = "talon-fat"

    resources:
        cpus = config["resources"]["bamCov"]["cpus"],
        time = config["resources"]["bamCov"]["time"],
        mem = config["resources"]["bamCov"]["mem"],


    shell:"""
        bamCoverage -p {resources.cpus} -b {input.bam} --normalizeUsing CPM -of bigwig -o {output.bw}
        """

