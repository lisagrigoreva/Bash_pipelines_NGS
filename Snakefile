import pandas as pd
# Import data_table
SAMPLES_INFO=pd.read_csv('data_table.tsv',sep="\t")
# Create column with sample path for the link
SAMPLES_INFO['Sample'] = SAMPLES_INFO["GSM"] + '_' + SAMPLES_INFO["Cell"] + '_' + SAMPLES_INFO["Target"]
SAMPLES_INFO['File'] = SAMPLES_INFO['File'].str.strip()
#Provide config file with genomes
configfile: "config.yaml"

# Target rule for all of the rules. Genome as wildcards
rule all:
    input:
       expand("qc/fastqc/{sample}.html", sample=SAMPLES_INFO['Sample']),
        expand("qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES_INFO['Sample']),
        "qc/multiqc/reads.html",
        expand("indexes/{genome}/{genome}.fa.gz", genome = config['genome']),
        expand("indexes/{genome}/{genome}.1.bt2", genome = config['genome']),
        expand("indexes/{genome}/{genome}.2.bt2",genome = config['genome']),
        expand("indexes/{genome}/{genome}.3.bt2", genome = config['genome']),
        expand("indexes/{genome}/{genome}.4.bt2",genome = config['genome']),
        expand("indexes/{genome}/{genome}.rev.1.bt2", genome = config['genome']),
        expand("indexes/{genome}/{genome}.rev.2.bt2", genome = config['genome']),
        expand("bams/{sample}_{genome}.bam", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        "qc/multiqc/bams.html",
        expand("bams_sorted/{sample}_{genome}.sorted.bam", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        expand("bams_sorted/{sample}_{genome}.sorted.bai", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        expand("bams_sorted/{sample}_{genome}.sorted.coverage.bw", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        "results.tar.gz"

#Fastqc rule that read on each file (with wrapper)
rule run_fastqc:
    input:
        lambda wildcards: expand("reads/{file_name}",
                                 file_name=SAMPLES_INFO.loc[SAMPLES_INFO['Sample'] == wildcards.sample, 'File'])
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
           "0.57.0/bio/fastqc"
#
# #Rule for multiqc with wrapper
rule multiqc:
    input:
        expand("qc/fastqc/{sample}.html", sample=SAMPLES_INFO['Sample'])
    output:
        "qc/multiqc/reads.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.57.0/bio/multiqc"
# Download genome
rule download:
    output:
        expand("indexes/{genome}/{genome}.fa.gz", genome = config['genome'])
    shell:
        "wget http://hgdownload.soe.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.fa.gz -P indexes/{config[genome]}/"

#Rule that build Bowtie2 indexes for genome
rule Bowtie2:
    input:
        rules.download.output
    output:
        "indexes/{genome}/{genome}.1.bt2",
        "indexes/{genome}/{genome}.2.bt2",
        "indexes/{genome}/{genome}.3.bt2",
        "indexes/{genome}/{genome}.4.bt2",
        "indexes/{genome}/{genome}.rev.1.bt2",
        "indexes/{genome}/{genome}.rev.2.bt2"
    conda:
        "envs/bowtie2.yml"
    params:
        output_name = expand("indexes/{genome}/{genome}", genome=config["genome"])
    shell:
            """
            gunzip -c {input} > {input}.tmp
          bowtie2-build {input}.tmp {params.output_name}
            """
# Alignment rule with wrapper
rule alignment:
    input:
        rules.Bowtie2.output,
        sample=lambda wildcards: expand("reads/{file_name}",
                                 file_name=SAMPLES_INFO.loc[SAMPLES_INFO['Sample'] == wildcards.sample, 'File'])
    output:
        "bams/{sample}_{genome}.bam"
    log:
        "logs/alignment/{sample}_{genome}.log"
    params:
        index = expand("indexes/{genome}/{genome}", genome=config["genome"]),
        extra=""
    threads: 4
    wrapper:
        "0.58.0/bio/bowtie2/align"

#Mutliqc with shell section
rule multiqc2:
    input:  expand(rules.alignment.log, sample=SAMPLES_INFO['Sample'], genome=config["genome"])
    log:
        "logs/multiqc2.log"
    output:
        "qc/multiqc/bams.html"
    conda:
        "envs/multiqc.yml"
    shell:
         "multiqc -d {input} -n {output}"
#Sort bam files using samtools
rule samtools_sort:
    input:
        expand("bams/{sample}_{genome}.bam", sample=SAMPLES_INFO['Sample'], genome=config["genome"])
    output:
        "bams_sorted/{sample}_{genome}.sorted.bam"
    log:
        "logs/samtools_sort/{sample}_{genome}.log"
    params:
        "-m 4G"
    threads:
        4
    wrapper:
        "0.59.2/bio/samtools/sort"
#Index bam files using samtools
rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "bams_sorted/{sample}_{genome}.sorted.bai"
    log:  "logs/samtools_index/{sample}_{genome}.log"
    wrapper:
        "0.59.2/bio/samtools/index"

#Coverage for each bam file
rule bigwig:
    input:
         bai = rules.samtools_index.output,
         bam = rules.samtools_sort.output
    output:
          "bams_sorted/{sample}_{genome}.sorted.coverage.bw"
    log: "logs/bigwig/{sample}_{genome}.log"
    conda:
          "envs/bamCoverage.yaml"
    shell:
         """
         bamCoverage \
             --bam {input.bam} \
             --outFileName {output} \
             --numberOfProcessors 2
                   """
 #Tar files
rule archive:
    input:
         fastqc = expand("qc/fastqc/{sample}.html", sample=SAMPLES_INFO['Sample']),
         zip = expand("qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES_INFO['Sample']),
         multiqc = "qc/multiqc/reads.html",
         bams = "qc/multiqc/bams.html",
         bc = expand("bams_sorted/{sample}_{genome}.sorted.coverage.bw", sample=SAMPLES_INFO['Sample'],genome=config["genome"])
    output:
          "results.tar.gz"
    log:
        "logs/tar.log"
    shell:
         "tar czvf {output} {input}"
























