import os.path


def flatten(input):
    flat_input = []
    for lst in input:
        if type(lst) is list:
            # If the element is of type list, iterate through the sublist
            for path in lst:
                if type(path) is list:
                    for sub in path:
                        flat_input.append(sub)
                else:
                    flat_input.append(path)         
        else:
            flat_input.append(lst)
    return (flat_input)


def get_fastqc_input(wildcards):
    if '_Trimmed' in wildcards.fname:
        fqc_input = "{outdir}/Trimmed_fastq/{{fname}}.{fqsuffix}"
    else:
        fqc_input = "{fq_dir}/{{fname}}.{fqsuffix}"

    return sorted(expand(fqc_input, **config))

rule FastQC:
    input: 
        get_fastqc_input
                  
    output:
         html = "{outdir}/FastQC/{{fname}}_fastqc.html".format(**config),
         zip = "{outdir}/FastQC/{{fname}}_fastqc.zip".format(**config)

    params:
        dir ="{outdir}/FastQC".format(**config),
        partition = "talon"

    log:
        "{outdir}/logs/fastqc/{{fname}}.log".format(**config)

    shell: "fastqc -o {params.dir} {input}"


# rule Trimmed_FastQC:
#     input: 
#         "{outdir}/Trimmed_fastq/{{fname}}.{fqsuffix}".format(**config)

                  
#     output:
#          html = "{outdir}/Trimmed_FastQC/{{fname}}_report.html".format(**config),
#          zip = "{outdir}/Trimmed_FastQC/{{fname}}_fastqc.zip".format(**config)

#     params:
#         dir = "{outdir}/Trimmed_FastQC".format(**config),
#         partition = "talon"

#     log:
#         "{outdir}/logs/Trimmedfastqc/{{fname}}.log".format(**config)

#     shell: "fastqc -o {params.dir} {input}" 


#rule combine_multiqc_files:
#    input: get_multiqc_input

def get_multiqc_input():
    input = []
    ## FastQC of Raw fastq files
    input.append([expand(f"{{outdir}}/FastQC/{sample}_{{fqext}}_fastqc.zip",**config) for sample in samples.index])

    ## FastQC of Trimmed fastq files
    if 'replicate' in samples:
            for replicate in set(samples['replicate']):
                input.append([expand(f"{{outdir}}/FastQC/{replicate}_{{fqext}}_Trimmed_fastqc.zip",**config)])
                input.append([expand(f"{{outdir}}/logs/Trimomatic/{replicate}_trimmomatic.out",**config)])

    else:
        input.append([expand(f"{{outdir}}/FastQC/{sample}_{{fqext}}_Trimmed_fastqc.zip",**config) for sample in samples.index])
        input.append([expand(f"{{outdir}}/logs/Trimomatic/{sample}_trimmomatic.out",**config) for sample in samples.index])

        

    if config["aligner"] == "STAR":
        # get alignment stats
        if 'replicate' in samples:
            for replicate in set(samples['replicate']):
                input.append([expand(f"{{outdir}}/Alignments/{replicate}_Log.final.out",**config)])
            
        else:
            input.append([expand(f"{{outdir}}/Alignments/{sample}_Log.final.out",**config) for sample in samples.index])
        

    elif config["aligner"] == "hisat2":
        if 'replicate' in samples:
            for replicate in set(samples['replicate']):
                input.append([expand(f"{{outdir}}/Alignments/{replicate}_align_stat.txt",**config)])
            
        else:
             input.append([expand(f"{{outdir}}/Alignments/{sample}_align_stat.txt",**config) for sample in samples.index])
        

    ### get featureCounts output if run
    stats = "{outdir}/Counts/geneCounts.txt.summary".format(**config)    
    if os.path.isfile(stats):
        input.append(stats)


    ## flatten list of inputs
    
    return(flatten(input))


rule combine_qc_files:
    input: get_multiqc_input()

    output:
        temp(expand("{outdir}/multiqc.tmp.files", **config)),

    params:
        partition = "talon"

    run:
        with open(output[0], mode="w") as out:
            out.write('\n'.join(input))

rule Multiqc:
    input:
        files=rules.combine_qc_files.output,

    output:
        "{outdir}/Multiqc_report.html".format(**config)   

    params:
        partition = "talon", 
        outdir = "{outdir}".format(**config)

    shell: """
        multiqc $(< {input.files}) -o {params.outdir} -n Multiqc_report.html
    """
