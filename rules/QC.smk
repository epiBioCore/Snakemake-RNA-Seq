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