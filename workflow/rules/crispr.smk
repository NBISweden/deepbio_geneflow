rule minced:
    input:
        opj("results", "assembly", "{assembly}", "contigs.fasta.gz")
    output:
        txt = opj("results", "minced", "{assembly}", "{assembly}.crisprs.txt.gz"),
        gff = opj("results", "minced", "{assembly}", "{assembly}.crisprs.gff.gz")
    params:
        tmp = opj("$TMPDIR", "{assembly}.minced"),
        fa = opj("$TMPDIR", "{assembly}.minced", "contigs.fasta"),
        txt = opj("$TMPDIR", "{assembly}.minced", "{assembly}.crisprs.txt"),
        gff = opj("$TMPDIR", "{assembly}.minced", "{assembly}.crisprs.gff"),
        outdir = lambda wildcards, output: os.path.dirname(output.txt)
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