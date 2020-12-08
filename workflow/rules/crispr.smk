rule minced:
    input:
        "results/assembly/{assembly}/final.contigs.fa"
    output:
        txt = report("results/minced/{assembly}/{assembly}.crisprs.txt.gz",
                     caption="../report/minced.rst", category="CRISPRs", subcategory="Minced"),
        gff = "results/minced/{assembly}/{assembly}.crisprs.gff.gz"
    params:
        tmp = "$TMPDIR/{assembly}.minced",
        fa = "$TMPDIR/{assembly}.minced/contigs.fasta",
        txt = "$TMPDIR/{assembly}.minced/{assembly}.crisprs.txt",
        gff = "$TMPDIR/{assembly}.minced/{assembly}.crisprs.gff",
        outdir = lambda wildcards, output: os.path.dirname(output.txt),
        account=config["project"]
    conda:
        "../envs/minced.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*2
    shell:
        """
        # Create tmpdir
        mkdir -p {params.tmp}
        # Unzip fasta file
        gunzip -c {input[0]} > {params.fa}
        # Run minced
        minced {params.fa} {params.txt} {params.gff}
        # Sync output and clean up
        gzip {params.txt} {params.gff}
        mv {params.tmp}/*.gz {params.outdir}
        rm -r {params.tmp}
        """