rule make_fasta_from_fastg:
    input:
        opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz")
    output:
        opj("results", "assembly", "{assembly}", "assembly_graph.nodes.fasta")
    params:
        tmp = opj("$TMPDIR", "{assembly}.fastg2fasta"),
        fastg = opj("$TMPDIR", "{assembly}.fastg2fasta", "assembly_graph.fastg")
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

rule bwa_index:
    input:
        opj("results", "assembly", "{assembly}", "assembly_graph.nodes.fasta")
    output:
        expand(opj("results", "assembly", "{{assembly}}",
                   "assembly_graph.nodes.fasta.{s}"),
               s = ["amb", "ann", "bwt", "pac", "sa"])
    log:
        opj("results", "logs", "plasmids", "{assembly}.bwa_index.log")
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input[0]} > {log} 2>&1
        """

rule bwa_mem:
    input:
        fasta = opj("results", "assembly", "{assembly}", "assembly_graph.nodes.fasta"),
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2"),
        index = expand(opj("results", "assembly", "{{assembly}}",
                   "assembly_graph.nodes.fasta.{s}"),
               s = ["amb", "ann", "bwt", "pac", "sa"])
    output:
        opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam"),
        opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam.bai")
    log:
        opj("results", "logs", "plasmids", "{assembly}.bwa.log")
    params:
        tmp = opj("$TMPDIR", "{assembly}.bwa"),
        R1 = opj("$TMPDIR", "{assembly}.bwa", "R1.fastq"),
        R2 = opj("$TMPDIR", "{assembly}.bwa", "R2.fastq"),
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        
        # Concatenate input
        gunzip -c {input.R1} > {params.R1}
        gunzip -c {input.R2} > {params.R2}
        
        # Map with bwa
        bwa mem -t {threads} {input.fasta} {params.R1} {params.R2} 2>{log} | \
            samtools view -buS - | samtools view -bF 0x0900 - | \
            samtools sort - > {params.tmp}/reads_pe_primary.sort.bam 
        # Index
        samtools index {params.tmp}/reads_pe_primary.sort.bam
        mv {params.tmp}/reads_pe_primary.sort.bam* {params.outdir}
        """

rule scapp:
    input:
        fastg = opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz"),
        bam = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam"),
        bai = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam.bai"),
        kmer = opj("results", "assembly", "{assembly}", "kmer")
    output:
        touch(opj("results", "scapp", "{assembly}", "{assembly}.confident_cycs.fasta"))
    log:
        opj("results", "logs", "plasmids", "{assembly}.scapp.log")
    conda:
        "../envs/scapp.yaml"
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir = opj("$TMPDIR", "{assembly}.scapp"),
        account = config["project"]
    shell:
        """
        set +e
        # Create tmpdir
        mkdir -p {params.tmpdir}
        
        # Unzip fastg
        gunzip -c {input.fastg} > {params.tmpdir}/{wildcards.assembly}.fastg
        
        # Get kmer size
        k=$(cat {input.kmer})
        
        # Run SCAPP
        scapp -p {threads} -g {params.tmpdir}/{wildcards.assembly}.fastg -k $k \
            -b {input.bam} -o {params.outdir} > {log} 2>&1
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            echo "NO PLASMIDS FOUND" > {output[0]}
        fi 
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

rule recycler:
    input:
        bam = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam"),
        bai = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam.bai"),
        kmer = opj("results", "assembly", "{assembly}", "kmer"),
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
        k=$(cat {input.kmer})
        recycle.py -g {params.graph} -k $k -b {input.bam} \
            -o {params.outdir} > {log} 2>&1
        """
