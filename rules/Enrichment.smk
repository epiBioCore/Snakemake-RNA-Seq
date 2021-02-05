rule GO_enrichment:
    input:
        sig = "{outdir}/DESeq2/{{contrast}}_DEtable_SIG_genes.csv".format(**config),
        species_anno = get_bioc_pkg_path,

    output:
       object = "{outdir}/Enrichment/{{contrast}}_{{GO}}_enrichment.rda".format(**config)

    params:
        species = get_bioc_species_pkg,
        out = "{outdir}/Enrichment".format(**config),
        partition = "talon",
        h = config["enrichment_plots"]["heigth"], 
        w = config["enrichment_plots"]["width"],
        u = config["enrichment_plots"]["units"]   
   
    wildcard_constraints:
        GO="GO_.+"
    log:
        "{outdir}/logs/Enrichment/{{contrast}}_{{GO}}.out".format(**config)
    benchmark:
        "{outdir}/benchmarks/Enrichment/{{contrast}}_{{GO}}.out".format(**config)
            
    resources:
        cpus = config["resources"]["Enrich"]["cpus"],
        time = config["resources"]["Enrich"]["time"],
        mem = config["resources"]["Enrich"]["mem"],

    script: "{workdir}/Scripts/GO.R".format(**config)    


rule GSEA_enrichment:
    input:
        de = "{outdir}/DESeq2/{{contrast}}_DEtable_ALL_genes.csv".format(**config),
        species_anno = get_bioc_pkg_path,

    output:
       object = "{outdir}/Enrichment/{{contrast}}_GSEA_enrichment.rda".format(**config)

    params:
        species = get_bioc_species_pkg,
        out = "{outdir}/Enrichment".format(**config),
        partition = "talon",
        h = config["enrichment_plots"]["heigth"], 
        w = config["enrichment_plots"]["width"],
        u = config["enrichment_plots"]["units"] ,  
   
  
    log:
        "{outdir}/logs/Enrichment/{{contrast}}_GSEA.out".format(**config)
    benchmark:
        "{outdir}/benchmarks/Enrichment/{{contrast}}_GSEA.out".format(**config)
            
    resources:
        cpus = config["resources"]["Enrich"]["cpus"],
        time = config["resources"]["Enrich"]["time"],
        mem = config["resources"]["Enrich"]["mem"],

    script: "{workdir}/Scripts/GSEA.R".format(**config)    

rule KEGG_module_enrichment:
    input:
        sig = "{outdir}/DESeq2/{{contrast}}_DEtable_SIG_genes.csv".format(**config),
        species_anno = get_bioc_pkg_path,


    output:
       object = "{outdir}/Enrichment/{{contrast}}_KEGG_Module_enrichment.rda".format(**config)

    params:
        species = config["annotation"]["species"],
        species_anno = get_bioc_species_pkg,
        out = "{outdir}/Enrichment".format(**config),
        partition = "talon",
        h = config["enrichment_plots"]["heigth"], 
        w = config["enrichment_plots"]["width"],
        u = config["enrichment_plots"]["units"]   
   
  
    log:
        "{outdir}/logs/Enrichment/{{contrast}}_KEGG_module.out".format(**config)
    benchmark:
        "{outdir}/benchmarks/Enrichment/{{contrast}}_KEGG_module.out".format(**config)
            
    resources:
        cpus = config["resources"]["Enrich"]["cpus"],
        time = config["resources"]["Enrich"]["time"],
        mem = config["resources"]["Enrich"]["mem"],

    script: "{workdir}/Scripts/KEGG_Module.R".format(**config)    

rule KEGG_enrichment:
    input:
        sig = "{outdir}/DESeq2/{{contrast}}_DEtable_SIG_genes.csv".format(**config),
        species_anno = get_bioc_pkg_path,


    output:
       object = "{outdir}/Enrichment/{{contrast}}_KEGG_enrichment.rda".format(**config)

    params:
        species = config["annotation"]["species"],
        species_anno = get_bioc_species_pkg,
        out = "{outdir}/Enrichment".format(**config),
        partition = "talon",
        h = config["enrichment_plots"]["heigth"], 
        w = config["enrichment_plots"]["width"],
        u = config["enrichment_plots"]["units"]   
   
  
    log:
        "{outdir}/logs/Enrichment/{{contrast}}_KEGG.out".format(**config)
    benchmark:
        "{outdir}/benchmarks/Enrichment/{{contrast}}_KEGG.out".format(**config)
            
    resources:
        cpus = config["resources"]["Enrich"]["cpus"],
        time = config["resources"]["Enrich"]["time"],
        mem = config["resources"]["Enrich"]["mem"],

    script: "{workdir}/Scripts/KEGG.R".format(**config)  

rule wikipathways_enrichment:
    input:
        sig = "{outdir}/DESeq2/{{contrast}}_DEtable_SIG_genes.csv".format(**config),
        species_anno = get_bioc_pkg_path,


    output:
       object = "{outdir}/Enrichment/{{contrast}}_wikipathways_enrichment.rda".format(**config)

    params:
        species = config["annotation"]["species"],
        species_anno = get_bioc_species_pkg,
        out = "{outdir}/Enrichment".format(**config),
        partition = "talon",
        h = config["enrichment_plots"]["heigth"], 
        w = config["enrichment_plots"]["width"],
        u = config["enrichment_plots"]["units"]   
   
  
    log:
        "{outdir}/logs/Enrichment/{{contrast}}_wikipathways.out".format(**config)
    benchmark:
        "{outdir}/benchmarks/Enrichment/{{contrast}}_wikipathways.out".format(**config)
            
    resources:
        cpus = config["resources"]["Enrich"]["cpus"],
        time = config["resources"]["Enrich"]["time"],
        mem = config["resources"]["Enrich"]["mem"],

    script: "{workdir}/Scripts/wikipathways.R".format(**config)    

rule spia:
    input:
        sig = "{outdir}/DESeq2/{{contrast}}_DEtable_SIG_genes.csv".format(**config),
        species_anno = get_bioc_pkg_path,


    output:
       object = "{outdir}/Enrichment/{{contrast}}_spia.txt".format(**config)

    params:
        species = config["annotation"]["species"],
        species_anno = get_bioc_species_pkg,
        fdr = config["DESeq2"]["fdr"],
        out = "{outdir}/Enrichment".format(**config),
        partition = "talon"
         
   
  
    log:
        "{outdir}/logs/Enrichment/{{contrast}}_spia.out".format(**config)
    benchmark:
        "{outdir}/benchmarks/Enrichment/{{contrast}}_spia.out".format(**config)
            
    resources:
        cpus = config["resources"]["Enrich"]["cpus"],
        time = config["resources"]["Enrich"]["time"],
        mem = config["resources"]["Enrich"]["mem"],

    script: "{workdir}/Scripts/SPIA.R".format(**config)   