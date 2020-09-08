from scripts.common import get_assembly_files

rule metaspades:
    input:
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2")
    output:
        fasta = expand(opj("results", "assembly", "{{assembly}}", "{f}.fasta.gz"),
               f = ["contigs", "scaffolds"]),
        fastg = opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz"),
        kmer = opj("results", "assembly", "{assembly}", "kmer")
    log:
        opj("results", "logs", "assembly", "{assembly}.spades.log")
    params:
        tmp=opj("$TMPDIR","{assembly}.metaspades"),
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        account=config["project"]
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*8
    conda:
        "../envs/metaspades.yaml"
    shell:
        """
        # Create directories
        mkdir -p {params.tmp}
        
        # Merge input
        gunzip -c {input.R1} > {params.tmp}/R1.fq
        gunzip -c {input.R2} > {params.tmp}/R2.fq
        
        # Run metaspades
        spades.py --meta \
            -t {threads} -1 {params.tmp}/R1.fq -2 {params.tmp}/R2.fq \
            -o {params.tmp}
        
        # Compress output
        gzip {params.tmp}/scaffolds.fasta {params.tmp}/contigs.fasta {params.tmp}/assembly_graph.fastg
        # Move output from temporary directory
        mv {params.tmp}/*.gz {params.output_dir}
        mv {params.tmp}/spades.log {params.tmp}/params.txt {params.output_dir}
        # Clean up
        rm -r {params.tmp}
        
        # Extract max kmer from spades log
        k=$(egrep -A 1 "^Assembly parameters:" {log} | grep "k:" | egrep -o "[0-9]+\]" | sed 's/]//g')
        echo "$k" > {output.kmer}
        """
