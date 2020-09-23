## Configuration file
import os
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['genome', 'metatxt', 'ncores', 'FASTQ', 'fqext1', 'fqext2', 'fqsuffix', 'output', 'run_trimming']
for k in configvars:
	if k not in config:
		config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
	if str is None:
		str = ''
	return str

#config['genome'] = sanitizefile(config['genome'])
config['metatxt'] = sanitizefile(config['metatxt'])

## Read metadata
if not os.path.isfile(config["metatxt"]):
  sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

import pandas as pd
samples = pd.read_csv(config["metatxt"], sep='\t')

if not set(['names','type']).issubset(samples.columns):
  sys.exit("Make sure 'names' and 'type' are columns in " + config["metatxt"])


## Sanitize provided input and output directories
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

outputdir = getpath(config["output"])
FASTQdir = getpath(config["FASTQ"])

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		outputdir + "MultiQC/multiqc_report.html",
		expand(outputdir + "methextract/{sample}_pe.bismark.cov.gz", sample = samples.names.values.tolist())

rule setup:
	input:
		outputdir + "Rout/pkginstall_state.txt",
		outputdir + "Rout/softwareversions.done"


## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

rule runbismark:
	input:
		expand(outputdir + "bismark/{sample}_pe.bam", sample = samples.names.values.tolist())


## Print the versions of all software packages
rule softwareversions:
	output:
		touch(outputdir + "Rout/softwareversions.done")
	conda:
		"envs/environment.yaml"
	shell:
		"echo -n 'cutadapt ' && cutadapt --version; "
		"fastqc --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##
## Generate Bismark index
rule bismarkindex:
	input:
		genome = config["genome"]
	output:
		config["genome"] + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		config["genome"] + "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
	log:
		outputdir + "logs/bismark_index.log"
	benchmark:
		outputdir + "benchmarks/bismark_index.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"bismark --version >> {log}; "
		"bismark_genome_preparation --verbose {input.genome}"

## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_{sample}.log"
	benchmark:
		outputdir + "benchmarks/fastqc_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## FastQC, trimmed reads
rule fastqctrimmed:
	input:
		fastq = outputdir + "FASTQtrimmed/{sample}.fq.gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_trimmed_{sample}.log"
	benchmark:
		outputdir + "benchmarks/fastqc_trimmed_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"



# The config.yaml files determines which steps should be performed
def multiqc_input(wildcards):
	input = []
	input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
	input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
	input.extend(expand(outputdir + "bismark/{sample}_pe.bam", sample = samples.names.values.tolist()))
	
	if config["run_trimming"]:
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()))
	return input

## Determine the input directories for MultiQC depending on the config file
def multiqc_params(wildcards):
	param = [outputdir + "FastQC",
	outputdir + "bismark"]
	if config["run_trimming"]:
		param.append(outputdir + "FASTQtrimmed")
	return param

## MultiQC
rule multiqc:
	input:
		multiqc_input
	output:
		outputdir + "MultiQC/multiqc_report.html"
	params:
		inputdirs = multiqc_params,
		MultiQCdir = outputdir + "MultiQC"
	log:
		outputdir + "logs/multiqc.log"
	benchmark:
		outputdir + "benchmarks/multiqc.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"
		

## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgalorePE:
	input:
		fastq1 = FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
		fastq2 = FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz",
		outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz"
	params:
		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
	log:
		outputdir + "logs/trimgalore_{sample}.log"
	benchmark:
		outputdir + "benchmarks/trimgalore_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore --clip_R1 5 --clip_R2 5 --three_prime_clip_R1 5 --three_prime_clip_R2 5 "
		"-o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt --paired {input.fastq1} {input.fastq2}"

## ------------------------------------------------------------------------------------ ##
## Bismark mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with bismark

rule bismarkPE:
	input:
		config["genome"] + "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
		config["genome"] + "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
		fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
		fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "bismark/{sample}_pe.bam"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/bismark_{sample}.log"
	benchmark:
		outputdir + "benchmarks/bismark_{sample}.txt"
	params:
		genome = config["genome"],
		bismarkdir = outputdir + "bismark",
		samplename = "{sample}"
	conda:
		"envs/environment.yaml"
	shell:
		"bismark --version >> {log}; "
		"bismark -p {threads} -B {params.samplename} -o {params.bismarkdir} --genome {params.genome} "
		"-1 {input.fastq1} -2 {input.fastq2}"
		

rule extractMethylation:
	input:
		outputdir + "bismark/{sample}_pe.bam"
	output:
		outputdir + "methextract/{sample}_pe.bismark.cov.gz"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/methextract_{sample}.log"
	benchmark:
		outputdir + "benchmarks/methextract_{sample}.txt"
	params:
		genome = config["genome"],
		methdir = outputdir + "methextract"
	conda:
		"envs/environment.yaml"
	shell:
		"bismark --version >> {log}; "
		"bismark_methylation_extractor -p --cytosine_report --comprehensive "
		" --genome_folder {params.genome} --no_overlap --multicore {threads} --buffer_size 2G "
		" -o {params.methdir} {input}"
	

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
