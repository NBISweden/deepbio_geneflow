##### download example files #####

rule download_synthetic:
    """
    Download pre-made synthetic metagenome from Zenodo
    """
    output:
        R1 = temp("data/synthetic_1.fastq.gz"),
        R2 = temp("data/synthetic_2.fastq.gz")
    log:
        "results/logs/testdata/synthetic.log"
    params:
        tar = "examples/data/synthetic.tar.gz",
        url = "https://zenodo.org/record/3737112/files/synthetic.tar.gz?download=1",
        outdir = lambda wildcards, output: os.path.dirname(output.R1)
    shell:
         """
         curl -L -s -o {params.tar} {params.url}
         tar -C {params.outdir} -xf {params.tar}
         rm {params.tar}
         """

rule generate_testdata:
    """
    Use seqtk to subsample the synthetic metagenome into examples
    """
    input:
        "data/synthetic_{i}.fastq.gz"
    output:
        "data/testdata/{sample}_{unit}_R{i}.fastq.gz"
    log:
        "results/logs/testdata/{sample}_{unit}_R{i}.log"
    conda:
        "../envs/examples.yaml"
    params:
        example_dataset_size=config["example_dataset_size"]
    shell:
         """
         seqtk sample -s {wildcards.unit} {input} \
            {params.example_dataset_size} | gzip -c > {output}
         """
