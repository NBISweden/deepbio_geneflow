from snakemake.utils import validate
import pandas as pd
import platform
import os
from os.path import join as opj
from scripts.common import parse_samples

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3:4.8.2"

##### load and validate config #####

if os.path.exists("config/config.yaml"):
    configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml", set_default=True)

##### load and validate samples #####
df = pd.read_csv(config["sample_list"], sep="\t")
validate(df, schema="../schemas/samples.schema.yaml")

##### parse samples #####
from scripts.common import parse_samples
samples, assemblies = parse_samples(df)

##### workflow settings #####

wildcard_constraints:
    unit="\d+",
    pair="R[12]",

rule link_samples:
    input:
        R1 = lambda wildcards: samples[wildcards.sample][wildcards.unit]["R1"],
        R2 = lambda wildcards: samples[wildcards.sample][wildcards.unit]["R2"]
    output:
        R1 = opj("data", "stage", "{sample}_{unit}_R1.fastq.gz"),
        R2 = opj("data", "stage", "{sample}_{unit}_R2.fastq.gz")
    params:
        abs_R1_in = lambda wildcards, input: os.path.abspath(input.R1),
        abs_R2_in = lambda wildcards, input: os.path.abspath(input.R2),
        abs_R1_out = lambda wildcards, output: os.path.abspath(output.R1),
        abs_R2_out = lambda wildcards, output: os.path.abspath(output.R2)
    shell:
        """
        ln -s {params.abs_R1_in} {params.abs_R1_out}
        ln -s {params.abs_R2_in} {params.abs_R2_out}
        """