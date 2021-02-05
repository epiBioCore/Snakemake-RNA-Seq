from snakemake.utils import validate
from snakemake.logging import logger
import pandas as pd
import os.path
import glob
import os
from pandas_schema import Column, Schema
from pandas_schema.validation import MatchesPatternValidation, IsDistinctValidation
from snakemake.exceptions import TerminatedException


## inspired by https://github.com/vanheeringen-lab/seq2science and 
# https://github.com/snakemake-workflows/rna-seq-star-deseq2

configfile: "config.yaml"
workdir: config['workdir']

samples = pd.read_table(config["sample_sheet"], dtype=str)
    #validate(samples, schema="schemas/samples.schema.yaml")
    
#samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index

#workdir: config["outdir"]

#samples = samples_df.index.get_level_values(0)
#samples = samples_df.index.get_level_values(1)
#samples = list(samples_df.unit.unique)


## sanitize samples sheet
samples.columns = samples.columns.str.strip()

assert all([col[0:7] not in ["Unnamed", ''] for col in samples]), \
    (f"\nEncountered unnamed column in {config['sample_sheet']}.\n" +
     f"Column names: {str(', '.join(samples.columns))}.\n")

# use pandasschema for checking if samples file is filed out correctly
allowed_pattern = r'[A-Za-z0-9_.\-%]+'
distinct_columns = ["sample"]
if "descriptive_name" in samples.columns:
    distinct_columns.append("descriptive_name")

distinct_schema = Schema(
    [Column(col, [MatchesPatternValidation(allowed_pattern),
                  IsDistinctValidation()] if col in distinct_columns else [MatchesPatternValidation(allowed_pattern)], allow_empty=True) for col in
     samples.columns])

errors = distinct_schema.validate(samples)

if len(errors):
    logger.error("\nThere are some issues with parsing the samples file:")
    for error in errors:
        logger.error(error)
    logger.error("")  # empty line
    raise TerminatedException

# for each column, if found in samples.tsv:
# 1) if it is incomplete, fill the blanks with replicate/sample names
# (sample names if replicates are not found/applicable)
# 2) drop column if it is identical to the replicate/sample column, or if not needed
if 'replicate' in samples:
    samples['replicate'] = samples['replicate'].mask(pd.isnull, samples['sample'])
    if samples['replicate'].tolist() == samples['sample'].tolist() or config.get('technical_replicates') == 'keep':
        samples = samples.drop(columns=['replicate'])
if 'condition' in samples:
    samples['condition'] = samples['condition'].mask(pd.isnull, samples['replicate']) if 'replicate' in samples else \
        samples['condition'].mask(pd.isnull, samples['sample'])
    if samples['condition'].tolist() == samples['sample'].tolist() or config.get('biological_replicates') == 'keep':
        samples = samples.drop(columns=['condition'])
if 'descriptive_name' in samples:
    samples['descriptive_name'] = samples['descriptive_name'].mask(pd.isnull, samples['replicate']) if \
        'replicate' in samples else samples['descriptive_name'].mask(pd.isnull, samples['sample'])
    if ('replicate' in samples and samples['descriptive_name'].to_list() == samples['replicate'].to_list()) or \
        samples['descriptive_name'].to_list() == samples['sample'].to_list():
        samples = samples.drop(columns=['descriptive_name'])
# if 'strandedness' in samples:
#     samples['strandedness'] = samples['strandedness'].mask(pd.isnull, 'nan')
#     if config.get('ignore_strandedness', True) or not any([field in list(samples['strandedness']) for field in ['yes', 'forward', 'reverse', 'no']]):
#         samples = samples.drop(columns=['strandedness'])
# if 'colors' in samples:
#     samples['colors'] = samples['colors'].mask(pd.isnull, '0,0,0')  # nan -> black
#     samples['colors'] = [color_parser(c) for c in samples['colors']]  # convert input to HSV color
#     if not config.get('create_trackhub', False):
#         samples = samples.drop(columns=['colors'])

# if 'replicate' in samples:
#     # check if replicate names are unique between assemblies
#     r = samples[['assembly', 'replicate']].drop_duplicates().set_index('replicate')
#     for replicate in r.index:
#         assert len(r[r.index == replicate]) == 1, \
#             ("\nReplicate names must be different between assemblies.\n" +
#              f"Replicate name '{replicate}' was found in assemblies {r[r.index == replicate]['assembly'].tolist()}.")

# check if sample, replicate and condition names are unique between the columns
for idx in samples.index:
    if "condition" in samples:
        assert idx not in samples["condition"].values, f"sample names, conditions, and replicates can not overlap. " \
                                                       f"Sample {idx} can not also occur as a condition"
    if "replicate" in samples:
        assert idx not in samples["replicate"].values, f"sample names, conditions, and replicates can not overlap. " \
                                                       f"Sample {idx} can not also occur as a replicate"

if "condition" in samples and "replicate" in samples:
    for cond in samples["condition"]:
        assert cond not in samples["replicate"].values, f"sample names, conditions, and replicates can not overlap. " \
                                                        f"Condition {cond} can not also occur as a replicate"

# validate samples file
sample_schemas = [f for f in os.listdir("schemas") ]
for schema in sample_schemas:
    validate(samples, schema=f"./schemas/{schema}")

sanitized_samples = copy.copy(samples)

samples = samples.set_index('sample')
samples.index = samples.index.map(str)

## 
config['fqext'] = [config['fqext1'], config['fqext2']]
assert sorted(config['fqext'])[0] == config['fqext1'], \
    ("\nThe paired-end filename suffixes must be lexicographically ordered!\n" +
     f"Example suffixes: fqext1: R1, fqext2: R2\n" +
     f"Your suffixes:    fqext1: {config['fqext1']}, fqext2: {config['fqext2']}\n")


#names = samples.descriptive_name.unique().tolist()
#replicates = samples.replicate.unique().tolist()
# def get_fastq(wildcards):
#     if config["paired"]:
#          expand(
#              [
#                  f"{{fq_dir}}/{wildcards.sample}_{fqext}.{{fqsuffix}}"
#                  for fqext in **config
#         return([r1,r2])
#     else:
#         r1 = "{fq_dir}/{{sample}}_{fqext1}.{fqsuffix}".format(**config),
#         return(r1)


### include rules
include: "rules/Alignment.smk"
include: "rules/Trimming.smk"
include: "rules/QC.smk"
include: "rules/Quantification.smk"
include: "rules/merge_replicates.smk"
include: "rules/DEG.smk"
include: "rules/Enrichment.smk"
include: "rules/Tracks.smk"


### helper functions

## flattens nested list of input paths
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

def get_raw_fastqc():
    input = []
    input.append(expand("{outdir}/FastQC/{{sample}}_{fqext1}.fastqc.zip".format(**config),sample = samples.index))
    input.append(expand("{outdir}/FastQC/{{sample}}_{fqext2}.fastqc.zip".format(**config),sample = samples.index))

    return(flatten(input))

def get_trimmed_fastqc():
    input = []    
    ## FastQC of Trimmed fastq files
    if 'replicate' in samples:
            for replicate in set(samples['replicate']):
                input.append([expand(f"{{outdir}}/Trimmed_FastQC/{replicate}_{{fqext}}_zip",**config)])

    else:
        input.append([expand(f"{{outdir}}/Trimmed_FastQC/{sample}_{{fqext}}_fastqc.zip",**config) for sample in samples.index])

    return(flatten(input))

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

def get_bigwigs():
    bws = set()
        

    # get all the bams
    if 'replicate' in samples:
        for replicate in set(samples['replicate']):
            bws.update([expand(f"{{outdir}}/Bigwigs/{replicate}.bw", **config)[0]])
            
    else:
        # suffix is defined in configuration_generic (noting or _custom)
        bws.update([expand(f"{{outdir}}/Bigwigs/{sample}.bw", **config)[0] for sample in samples.index])
        
    return bws


## fix rule all
rule all:
    input:
         expand("{outdir}/Enrichment/{contrast}_{GO}_enrichment.rda",contrast=get_contrast(),GO=["GO_BP","GO_CC","GO_MF"],outdir=config["outdir"]),
         expand("{outdir}/Enrichment/{contrast}_GSEA_enrichment.rda",contrast=get_contrast(),outdir=config["outdir"]),
         expand("{outdir}/Enrichment/{contrast}_KEGG_enrichment.rda",contrast=get_contrast(),outdir=config["outdir"]),
         expand("{outdir}/Enrichment/{contrast}_KEGG_Module_enrichment.rda",contrast=get_contrast(),outdir=config["outdir"]),
         expand("{outdir}/Enrichment/{contrast}_wikipathways_enrichment.rda",contrast=get_contrast(),outdir=config["outdir"]),
         expand("{outdir}/Enrichment/{contrast}_spia.txt",contrast=get_contrast(),outdir=config["outdir"]),
         expand("{outdir}/DESeq2/PCA_{dim}_colored_by_{var}.pdf",dim=["PC1vsPC2","PC1vsPC3","PC2vsPC3"],var=config["pca"]["col"],outdir=config["outdir"]),
         expand("{outdir}/DESeq2/dist_clustering.pdf",contrast=get_contrast(),outdir=config["outdir"]),
         expand("{outdir}/DESeq2/{contrast}_DEtable_ALL_genes.csv",contrast=get_contrast(),outdir=config["outdir"]),
         "{outdir}/DESeq2/DE_summary.csv".format(outdir=config["outdir"]),
         "{outdir}/Multiqc_report.html".format(outdir=config["outdir"]),
         get_bigwigs()


rule align_all:
    """
    align each sample against its assembly
    """
    input:
         get_bams()['map'],
         get_bams()['ind'],
         "{outdir}/Multiqc_report.html".format(outdir=config["outdir"]),


rule fpkms_all:
    input: 
        "{outdir}/Cuffnorm/genes.fpkm_table".format(**config),
        "{outdir}/Multiqc_report.html".format(outdir=config["outdir"]),


rule counts_all:
    input:
        "{outdir}/Counts/geneCounts.txt".format(**config),
        "{outdir}/Multiqc_report.html".format(outdir=config["outdir"])


        
rule bigwigs_all:
    input:
        get_bigwigs(),
        "{outdir}/Multiqc_report.html".format(outdir=config["outdir"]),

rule DESeq2_all:
    input:
        get_bigwigs(),
        "{outdir}/Multiqc_report.html".format(outdir=config["outdir"]),
        "{outdir}/DESeq2/DE_summary.csv".format(outdir=config["outdir"]),
        expand("{outdir}/DESeq2/PCA_{dim}_colored_by_{var}.pdf",dim=["PC1vsPC2","PC1vsPC3","PC2vsPC3"],var=config["pca"]["col"],outdir=config["outdir"]),
        expand("{outdir}/DESeq2/dist_clustering.pdf",contrast=get_contrast(),outdir=config["outdir"]),
        expand("{outdir}/DESeq2/{contrast}_DEtable_ALL_genes.csv",contrast=get_contrast(),outdir=config["outdir"]),

rule PCA_all:
    input:
        expand("{outdir}/DESeq2/PCA_{dim}_colored_by_{var}.pdf",dim=["PC1vsPC2","PC1vsPC3","PC2vsPC3"],var=config["pca"]["col"],outdir=config["outdir"]),

