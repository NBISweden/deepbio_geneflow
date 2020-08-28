rule scapp:
    input:
        graph = opj("results", "assembly", "{assembly}", "assembly_graph.fastg"),
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2"),
        log = opj("results", "logs", "assembly", "{assembly}.spades.log")
    output:
        touch(opj("results", "scapp", "{assembly}", "scapp.done"))
    conda:
        "../envs/scapp.yaml"
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir = opj("$TMPDIR", "scapp.{assembly}")
    shell:
        """
        # Create tmpdir
        mkdir -p {params.tmpdir}
        
        # Extract max kmer from spades log
        egrep -A 1 "^Assembly parameters:" {input.log} | grep "k:" | egrep -o "[0-9]+\]" | sed 's/]//g'
         
        # Concatenate input reads
        gunzip -c {input.R1} > {params.tmpdir}/R1.fq
        gunzip -c {input.R2} > {params.tmpdir}/R2.fq
        
        # Run SCAPP
        scapp -p {threads} -g {input.graph} \
            -r1 {params.tmpdir}/R1.fq -r2 {params.tmpdir}/R2.fq \
            -o {params.outdir}
        """
