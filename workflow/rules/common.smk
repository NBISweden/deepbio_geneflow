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
        lambda wildcards: samples[wildcards.sample][wildcards.unit][wildcards.R]
    output:
        opj("results", "stage", "{sample}_{unit}_{R}.fastq.gz")
    params:
        abs_in = lambda wildcards, input: os.path.abspath(input[0]),
        abs_out = lambda wildcards, output: os.path.abspath(output[0]),
    shell:
        """
        ln -s {params.abs_in} {params.abs_out}
        """