rule scapp:
    input:
        graph = opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz"),
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2"),
        log = opj("results", "logs", "assembly", "{assembly}.spades.log")
    output:
        touch(opj("results", "scapp", "{assembly}", "scapp.done"))
    log:
        opj("results", "logs", "plasmids", "{assembly}.scapp.log")
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
        k=$(egrep -A 1 "^Assembly parameters:" {input.log} | grep "k:" | egrep -o "[0-9]+\]" | sed 's/]//g')
         
        # Concatenate input reads
        gunzip -c {input.R1} > {params.tmpdir}/R1.fq
        gunzip -c {input.R2} > {params.tmpdir}/R2.fq
        
        # Unzip graph
        gunzip -c {input.graph} > {params.tmpdir}/graph.fastg
        
        # Run SCAPP
        scapp -p {threads} -g {params.tmpdir}/graph.fastg -k $k \
            -r1 {params.tmpdir}/R1.fq -r2 {params.tmpdir}/R2.fq \
            -o {params.outdir} > {log} 2>&1
        """
