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
        tmpdir = opj("$TMPDIR", "scapp.{assembly}"),
        account = config["project"]
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

rule mpspades:
    input:
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2")
    output:
        expand(opj("results", "mpspades", "{{assembly}}", "{f}.fasta.gz"),
               f = ["contigs", "scaffolds"]),
        touch(opj("results", "mpspades", "{assembly}", "assembly_graph.fastg.gz"))
    log:
        opj("results", "logs", "plasmids", "{assembly}.mpspades.log")
    threads: 4
    params:
        tmp=opj("$TMPDIR","{assembly}.mpspades"),
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        account=config["project"]
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
        spades.py --meta --plasmid \
            -t {threads} -1 {params.tmp}/R1.fq -2 {params.tmp}/R2.fq \
            -o {params.tmp} > {log} 2>&1
        
        # Move output
        # Compress output
        gzip {params.tmp}/scaffolds.fasta {params.tmp}/contigs.fasta
        # Move output from temporary directory
        mv {params.tmp}/*.gz {params.output_dir}
        mv {params.tmp}/spades.log {params.tmp}/params.txt {params.output_dir}
        # Clean up
        rm -r {params.tmp}
        """

rule make_fasta_from_fastg:
    input:
        opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz")
    output:
        opj("results", "recycler", "{assembly}", "assembly_graph.nodes.fasta")
    params:
        tmp = opj("$TMPDIR", "{assembly}.recycler"),
        fastg = opj("$TMPDIR", "{assembly}.recycler", "assembly_graph.fastg")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60
    conda:
        "../envs/recycler.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        gunzip -c {input[0]} > {params.fastg}
        make_fasta_from_fastg.py -g {params.fastg} -o {output[0]}
        rm -r {params.tmp}
        """

rule recycler_index:
    input:
        opj("results", "recycler", "{assembly}", "assembly_graph.nodes.fasta")
    output:
        expand(opj("results", "recycler", "{{assembly}}",
                   "assembly_graph.nodes.fasta.{s}"),
               s = ["amb", "ann", "bwt", "pac", "sa"])
    log:
        opj("results", "logs", "plasmids", "{assembly}.bwa_index.log")
    conda:
        "../envs/recycler.yaml"
    shell:
        """
        bwa index {input[0]} > {log} 2>&1
        """

rule recycler_map:
    input:
        fasta = opj("results", "recycler", "{assembly}", "assembly_graph.nodes.fasta"),
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2"),
        index = expand(opj("results", "recycler", "{{assembly}}",
                   "assembly_graph.nodes.fasta.{s}"),
               s = ["amb", "ann", "bwt", "pac", "sa"])
    output:
        opj("results", "recycler", "{assembly}", "reads_pe_primary.sort.bam"),
        opj("results", "recycler", "{assembly}", "reads_pe_primary.sort.bam.bai")
    log:
        opj("results", "logs", "plasmids", "{assembly}.bwa.log")
    params:
        tmp = opj("$TMPDIR", "{assembly}.recycler"),
        R1 = opj("$TMPDIR", "{assembly}.recycler", "R1.fastq"),
        R2 = opj("$TMPDIR", "{assembly}.recycler", "R2.fastq"),
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/recycler.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        
        # Concatenate input
        gunzip -c {input.R1} > {params.R1}
        gunzip -c {input.R2} > {params.R2}
        
        # Map with bwa
        bwa mem -t {threads} {input.fasta} {params.R1} {params.R2} 2>{log} | \
            samtools view -buS - | samtools view -bF 0x0800 - | \
            samtools sort - > {params.tmp}/reads_pe_primary.sort.bam 
        # Index
        samtools index {params.tmp}/reads_pe_primary.sort.bam
        mv {params.tmp}/reads_pe_primary.sort.bam* {params.outdir}
        """

rule recycler:
    input:
        bam = opj("results", "recycler", "{assembly}", "reads_pe_primary.sort.bam"),
        bai = opj("results", "recycler", "{assembly}", "reads_pe_primary.sort.bam.bai"),
        log = opj("results", "logs", "assembly", "{assembly}.spades.log"),
        graph = opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz")
    output:
        touch(opj("results", "recycler", "{assembly}", "done"))
    log:
        opj("results", "logs", "plasmids", "{assembly}.recycler.log")
    params:
        tmp = opj("$TMPDIR", "{assembly}.recycler"),
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        graph = opj("$TMPDIR", "{assembly}.recycler", "{assembly}.fastg")
    conda:
        "../envs/recycler.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        gunzip -c {input.graph} > {params.graph}
        # Extract max kmer from spades log
        k=$(egrep -A 1 "^Assembly parameters:" {input.log} | grep "k:" | egrep -o "[0-9]+\]" | sed 's/]//g')
        recycle.py -g {params.graph} -k $k -b {input.bam} \
            -o {params.outdir} > {log} 2>&1
        """
