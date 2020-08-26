localrules: generate_metaspades_input

rule generate_metaspades_input:
    """Generate input files for use with Metaspades"""
    input:
        lambda wildcards: get_assembly_files(assemblies[wildcards.assembly])
    output:
        R1=temp(opj("results", "assembly", "{assembly}","R1.fq")),
        R2=temp(opj("results", "assembly", "{assembly}","R2.fq"))
    params:
        assembly = lambda wildcards: assemblies[wildcards.assembly]
    script:
        "../scripts/assembly_utils.py"

rule metaspades:
    input:
        R1=temp(opj("results", "assembly", "{assembly}","R1.fq")),
        R2=temp(opj("results", "assembly", "{assembly}","R2.fq"))
    output:
        opj("results", "assembly", "{assembly}", "final_contigs.fa")
    log:
        opj("results", "logs", "assembly", "{assembly}.spades.log")
    params:
        tmp=opj("$TMPDIR","{assembly}.metaspades"),
        output_dir=opj("results", "assembly", "{assembly}")
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/metaspades.yml"
    shell:
        """
        # Create directories
        mkdir -p {params.tmp}
        
        # Run metaspades
        metaspades.py \
            -t {threads} -1 {input.R1} -2 {input.R2} \
            -o {params.tmp} > {log} 2>&1
        
        # Move output from temporary directory        
        mv {params.tmp}/* {params.output_dir}
        mv {params.output_dir}/scaffolds.fasta {params.output_dir}/final_contigs.fa       
        """