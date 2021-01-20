




def get_contrast():
    contrasts = []

    c = pd.read_table(config["comparisons_sheet"],delimiter = "\t")
    ## comparisons file should have 4 columns:
    ## Batch/covariate Column Group1 Group2
    ## where Column is the column name in sample sheet containing groups 1 and 2 
    ## or 3 columns Column Group1 Group2
    if len(c.columns) == 4:
        for index,row in c.iterrows():
            batch = row[0].strip()
            col = row[1].strip()
            g1 = row[2].strip()
            g2 = row[3].strip()
            ## if contrasts with batch and without batch are in same sheet
            if batch == "":
                contrast = "%s_%s_vs_%s" % (col,g1,g2)
            else:
                 contrast = "%s+%s_%s_vs_%s" % (batch,col,g1,g2)

            contrasts.append(contrast)        

    elif len(c.columns) == 3:
        for index,row in c.iterrows():
            col = row[0].strip()
            g1 = row[1].strip()
            g2 = row[2].strip()
            contrast = "%s_%s_vs_%s" % (col,g1,g2)

            contrasts.append(contrast) 

    return(contrasts)               


def get_bioc_species_pkg(wildcards):
     """Get the package bioconductor package name for the the species in config.yaml"""
     species_letters = config["annotation"]["species"][0:2].capitalize()
     return "org.{species}.eg.db".format(species=species_letters)

# def get_bioc_TxDb_pkg(wildcards):
#     """Get the package bioconductor package name for the the species in config.yaml"""
#     species = config["txdb_annotation"]["species"].capitalize()
#     annotation = config["txdb_annotation"]["annotation"]
#     build = config["txdb_annotation"]["build"]
#     version = config["txdb_annotation"]["build"]

#     if config["txdb_annotation"]["annotation"] == "UCSC":
#         pkg = "TxDb.{species}.{annotation}.{build}.knownGene".format(species = species,annotation = annotation,build=build)
#     elif config["txdb_annotation"]["annotation"] == "Ensembl":
#         pkg = "EnsDb.{species}.v{version}".format(species=species,version=version) 
#     return pkg       

def get_bioc_pkg_path(wildcards):
    return "resources/bioconductor_anno/lib/R/library/{pkg}".format(pkg=get_bioc_species_pkg(wildcards))



localrules: all, download_bioconductor_annotation_packages





rule DESeq2_init:
    input:
        counts = "{outdir}/Counts/geneCounts_for_DESEq2.csv".format(**config),
        cuff_fpkms = "{outdir}/Cuffnorm/genes.fpkm_table".format(**config),
        len = "{outdir}/Counts/gene_lengths.csv".format(**config) ## check spelling

    output:
        dds = "{outdir}/DESeq2/DESeqDataSet_init.rda".format(**config),
        counts = "{outdir}/DESeq2/normed_counts.csv".format(**config),
        fpkms = "{outdir}/DESeq2/fpkm_values.csv".format(**config),


    benchmark:
        "{outdir}/benchmarks/DESeq2/deseq2_init.out".format(**config)    

    log:
        "{outdir}/logs/DESeq2/deseq2_init.out".format(**config)    

    params:
        meta = config["sample_sheet"],
        comparisons = config["comparisons_sheet"],
        replicates=True if "replicate" in samples else False,
        partition="talon"

    resources:
        cpus = config["resources"]["DESeq2"]["cpus"],
        time = config["resources"]["DESeq2"]["time"],
        mem = config["resources"]["DESeq2"]["mem"],

    script: "{workdir}/Scripts/DESeq2_init.R".format(**config)    

rule DESeq2:
    input:
        dds = "{outdir}/DESeq2/DESeqDataSet_init.rda".format(**config),
        counts = "{outdir}/DESeq2/normed_counts.csv".format(**config),
        fpkms = "{outdir}/DESeq2/fpkm_values.csv".format(**config),
        species_anno = get_bioc_pkg_path,

    output:
        de = "{outdir}/DESeq2/{{contrast}}_DEtable_ALL_genes.csv".format(**config),
        ma = "{outdir}/DESeq2/{{contrast}}_maPlot.pdf".format(**config)
    params:
        org = config["annotation"],
        fdr = config["DESeq2"]["fdr"],
        species = get_bioc_species_pkg,
        ids = config["annotation"]["ids"],
        partition = "talon"

    log:
        "{outdir}/logs/DESeq2/DESeq2_{{contrast}}".format(**config)
    benchmark:
        "{outdir}/benchmarks/DESeq2/DESeq2_{{contrast}}".format(**config)
            
    resources:
        cpus = config["resources"]["DESeq2"]["cpus"],
        time = config["resources"]["DESeq2"]["time"],
        mem = config["resources"]["DESeq2"]["mem"],

    script: "{workdir}/Scripts/DESeq2.R".format(**config)    


rule download_bioconductor_annotation_packages:
    output:
        directory("resources/bioconductor_anno/lib/R/library/{package}")

    params:
        path=lambda wc, output: Path(output[0]).parents[3],
    
    log:
        "results/logs/R/install_{package}.txt"
    shell: """
     conda create --quiet --yes -p {params.path} --channel bioconda bioconductor-{wildcards.package}
     """

rule pca:
    input:
        dds = "{outdir}/DESeq2/DESeqDataSet_init.rda".format(**config),

    output:
        expand("{outdir}/DESeq2/PCA_{dim}_colored_by_{var}.pdf",dim=["PC1vsPC2","PC1vsPC3","PC2vsPC3"],var=config["pca"]["col"],outdir=config["outdir"]),
        rld = "{outdir}/DESeq2/rld.rda".format(**config)

    params:
        labs = config["pca"]["labels"],
        col = config["pca"]["col"],
        lab_size = config["pca"]["label_size"],
        partition="talon"

    log:
        "{outdir}/logs/DESeq2/PCA.log".format(**config)
    benchmark:
        "{outdir}/benchmarks/DESeq2/pca.benchmark".format(**config)
     

    script: "{workdir}/Scripts/PCA.R".format(**config)    

rule sample_clustering:
    input:
        dds = "{outdir}/DESeq2/DESeqDataSet_init.rda".format(**config)

    output:
        hm = "{outdir}/DESeq2/dist_clustering.pdf".format(**config)

    params:
        row_labels = config["sample_clustering"]["row_labels"],
        col_labels = config["sample_clustering"]["col_labels"],
        partition="talon"

    log:
        "{outdir}/logs/DESeq2/sampleHeatmap.log".format(**config)

    benchmark:
        "{outdir}/benchmarks/DESeq2/sampleHeatmap.benchmark".format(**config)
     

    script: "{workdir}/Scripts/Clustering.R".format(**config)    

rule DE_summary:
    input:
        expand("{outdir}/DESeq2/{{contrast}}_DEtable_ALL_genes.csv".format(**config),contrast=get_contrast())

    output:
        s = "{outdir}/DESeq2/DE_summary.csv".format(**config)

    params:
        fdr = config["DESeq2"]["fdr"],
        partition = "talon" 

    log:
        "{outdir}/logs/DESeq2/DE_summary.log".format(**config)

    benchmark:
        "{outdir}/benchmarks/DESeq2/DE_summary.benchmark".format(**config)

    script: "{workdir}/Scripts/DE_summary.R".format(**config)    
  