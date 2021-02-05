def get_bams():
    bamfiles = {}
    
    map_ext = 'bam'
    map_idx = 'bai'
    
    bamfiles['map'] = set(); bamfiles['ind'] = set()

    # get all the bams
    if 'replicate' in samples:
        for replicate in set(samples['replicate']):
            bamfiles['map'].update([expand(f"{{outdir}}/Alignments/{replicate}_sorted.{map_ext}", **config)[0]])
            bamfiles['ind'].update([expand(f"{{outdir}}/Alignments/{replicate}_sorted.{map_ext}.{map_idx}", **config)[0]])
            
    else:
        # suffix is defined in configuration_generic (noting or _custom)
        bamfiles['map'].update([expand(f"{{outdir}}/Alignments/{sample}_sorted.{map_ext}", **config)[0] for sample in samples.index])
        bamfiles['ind'].update([expand(f"{{outdir}}/Alignments/{sample}_sorted.{map_ext}.{map_idx}", **config)[0] for sample in samples.index])
        

    return bamfiles




rule cuffnorm:
    input: get_bams()['map'],

    output:
        dir = "{outdir}/Cuffnorm/genes.fpkm_table".format(**config),


    params:
        dir = directory("{outdir}/Cuffnorm".format(**config)),
        gtf= config["gtf"], 
        lib= config["cuffnorm"]["lib"],
        partition = "talon-fat"

    resources:
        cpus = config["resources"]["cuffnorm"]["cpus"],
        time = config["resources"]["cuffnorm"]["time"],
        mem = config["resources"]["cuffnorm"]["mem"],

    benchmark: 
        "{outdir}/benchmarks/cuffnorm.bench.txt".format(**config)
    
    log: 
        "{outdir}/logs/cuffnorm.out".format(**config)

    shell:"""
        cuffnorm -o {params.dir} -p {resources.cpus} --library-type {params.lib} {params.gtf} {input} 2> {log}
        """

rule featureCounts:
    input: get_bams()['map'],

    output:
        counts = "{outdir}/Counts/geneCounts.txt".format(**config),
        stats = "{outdir}/Counts/geneCounts.txt.summary".format(**config)    
    
    log:
       "{outdir}/logs/featureCounts.out".format(**config) 

    benchmark:
       "{outdir}/benchmarks/featureCounts.bench.txt".format(**config) 
   

    params:
        gtf = config["gtf"],
        strand = config["featureCounts"]["strand"],
        other = config["featureCounts"]["other"],
        partition = "talon-fat"

    resources:
        cpus = config["resources"]["featureCounts"]["cpus"],
        time = config["resources"]["featureCounts"]["time"],
        mem = config["resources"]["featureCounts"]["mem"],


    shell: """
        featureCounts -T {resources.cpus} -a {params.gtf} -s {params.strand} {params.other} -o {output.counts} {input} 2> {log}
        """

rule prepfeatureCounts:
    input:
        counts = "{outdir}/Counts/geneCounts.txt".format(**config)

    output: 
        counts = "{outdir}/Counts/geneCounts_for_DESEq2.csv".format(**config),
        annot = "{outdir}/Counts/gene_lengths.csv".format(**config),
    
    params:
        partition = "talon"

    log:
        "{outdir}/logs/R/prepfeatureCounts.txt".format(**config)    

    script: "{workdir}/Scripts/prepfeatureCounts.R".format(**config)    

      