from scripts.common import get_assembly_files

rule generate_metaspades_input:
    """Generate input files for use with Metaspades"""
    input:
        lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], wildcards.R)
    output:
        fq = temp(opj("results", "assembly", "{assembly}","{R}.fq"))
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    params:
        assembly = lambda wildcards: assemblies[wildcards.assembly],
        R = lambda wildcards: wildcards.R,
        tmp_out = opj("$TMPDIR", "{assembly}.{R}.fq")
    script:
        "../scripts/assembly_utils.py"

rule metaspades:
    input:
        R1 = opj("results", "assembly", "{assembly}","R1.fq"),
        R2 = opj("results", "assembly", "{assembly}","R2.fq")
    output:
        expand(opj("results", "assembly", "{{assembly}}", "{f}.fasta"),
               f = ["contigs", "scaffolds"]),
        opj("results", "assembly", "{assembly}", "assembly_graph.fastg")
    log:
        opj("results", "logs", "assembly", "{assembly}.spades.log")
    params:
        tmp=opj("$TMPDIR","{assembly}.metaspades"),
        output_dir=lambda wildcards, output: os.path.dirname(output[0])
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/metaspades.yaml"
    shell:
        """
        # Create directories
        mkdir -p {params.tmp}
        
        # Run metaspades
        spades.py --meta \
            -t {threads} -1 {input.R1} -2 {input.R2} \
            -o {params.tmp} > {log} 2>&1
        
        # Move output from temporary directory        
        mv {params.tmp}/* {params.output_dir}       
        """
