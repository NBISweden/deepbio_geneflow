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
    pair="R[12]"

localrules: compress_testdata

rule compress_testdata:
    input:
        "data/testdata/test_R{i}.fastq"
    output:
        "data/testdata/test_R{i}.fastq.gz"
    shell:
        """
        gzip -c {input} > {output}
        """
