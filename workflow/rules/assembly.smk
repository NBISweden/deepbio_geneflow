from scripts.common import get_assembly_files

localrules:
    fasta2fastg,
    symlink_fasta

rule megahit:
    input:
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2")
    output:
        "results/assembly/{assembly}/final.contigs.fa"
    log:
        "results/assembly/{assembly}/log"
    params:
        account=config["project"],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        k = ",".join([str(k) for k in config["megahit"]["k"]]),
        R1 = lambda wildcards, input: ",".join([f for f in input.R1]),
        R2 = lambda wildcards, input: ",".join([f for f in input.R2]),
    conda:
        "../envs/megahit.yaml"
    threads: 20
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*48
    shell:
        """
        rm -rf {params.outdir}
        # Run megahit
        megahit -1 {params.R1} -2 {params.R2} -o {params.outdir} \
            --k-list {params.k} -t {threads}
        """

rule fasta2fastg:
    input:
        "results/assembly/{assembly}/final.contigs.fa"
    output:
        "results/assembly/{assembly}/final.contigs.fastg"
    params:
        indir = lambda wildcards, input: os.path.dirname(input[0])
    conda: "../envs/megahit.yaml"
    shell:
        """
        # Get k-max from log
        files=$(ls {params.indir}/intermediate_contigs/*.final.contigs.fa)
        k=$(basename -a $files | cut -f1 -d '.' | sed 's/k//g' | sort -n | tail -n 1)
        # Convert to fastg format
        megahit_toolkit contig2fastg $k {params.indir}/intermediate_contigs/k$k.contigs.fa > {output[0]}
        """

rule symlink_fasta:
    input:
        "results/assembly/{assembly}/final.contigs.fa"
    output:
        "results/assembly/{assembly}/assembly_graph.nodes.fasta"
    params:
        indir=lambda wildcards, input: os.path.dirname(input[0])
    shell:
        """
        # Get k-max from log
        files=$(ls {params.indir}/intermediate_contigs/*.final.contigs.fa)
        k=$(basename -a $files | cut -f1 -d '.' | sed 's/k//g' | sort -n | tail -n 1)
        # Symlink corresponding fasta file
        ln -s $(pwd)/{params.indir}/intermediate_contigs/k$k.contigs.fa {output[0]}
        """
