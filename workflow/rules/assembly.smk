from scripts.common import get_assembly_files

rule metaspades:
    input:
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2")
    output:
        expand(opj("results", "assembly", "{{assembly}}", "{f}.fasta"),
               f = ["contigs", "scaffolds"]),
        opj("results", "assembly", "{assembly}", "assembly_graph.fastg")
    log:
        opj("results", "logs", "assembly", "{assembly}.spades.log")
    params:
        tmp=opj("$TMPDIR","{assembly}.metaspades"),
        output_dir=lambda wildcards, output: os.path.dirname(output[0])
    threads: 20
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4,
        mem_mb=lambda wildcards, attempt: attempt*128000
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
            -o {params.tmp} > {log} 2>&1
        
        # Clean up input files
        rm {params.tmp}/R1.fq {params.tmp}/R2.fq
        # Move output from temporary directory
        mv {params.tmp}/* {params.output_dir}       
        """
