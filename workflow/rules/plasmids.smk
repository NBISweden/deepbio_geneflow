rule bowtie_index:
    input:
        "results/assembly/{assembly}/assembly_graph.nodes.fasta"
    output:
        expand("results/assembly/{{assembly}}/assembly_graph.nodes.{s}.bt2l",
               s = ["1","2","3","4","rev.1","rev.2"])
    log:
        "results/logs/plasmids/{assembly}.bowtie_index.log"
    params:
        account = config["project"],
        prefix = "results/assembly/{assembly}/assembly_graph.nodes"
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*2
    conda:
        "../envs/bowtie.yaml"
    threads: 10
    shell:
        """
        bowtie2-build -t {threads} --large-index {input} {params.prefix}
        """

rule bowtie2:
    input:
        fasta = "results/assembly/{assembly}/assembly_graph.nodes.fasta",
        R1 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R1"),
        R2 = lambda wildcards: get_assembly_files(assemblies[wildcards.assembly], "R2"),
        index = expand("results/assembly/{{assembly}}/assembly_graph.nodes.{s}.bt2l",
               s = ["1","2","3","4","rev.1","rev.2"])
    output:
        "results/assembly/{assembly}/reads_pe_primary.sort.bam",
        "results/assembly/{assembly}/reads_pe_primary.sort.bam.bai"
    log:
        "results/logs/plasmids/{assembly}.bwa.log"
    params:
        account=config["project"],
        tmp = "$TMPDIR/{assembly}.bowtie",
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        prefix = "results/assembly/{assembly}/assembly_graph.nodes"
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    conda:
        "../envs/bowtie.yaml"
    shell:
        """
        mkdir -p {params.tmp}
                
        # Map with bwa
        bowtie2 --very-sensitive -x {params.prefix} -p {threads} -1 {input.R1} -2 {input.R2} 2>{log} | \
            samtools view -buS - | samtools view -bF 0x0900 - | \
            samtools sort - > {params.tmp}/reads_pe_primary.sort.bam 
        # Index
        samtools index {params.tmp}/reads_pe_primary.sort.bam
        mv {params.tmp}/reads_pe_primary.sort.bam* {params.outdir}
        # Clean up
        rm -r {params.tmp}
        """

rule scapp:
    input:
        fastg = opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz"),
        bam = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam"),
        bai = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam.bai"),
        kmer = opj("results", "assembly", "{assembly}", "kmer")
    output:
        report(touch(opj("results", "scapp", "{assembly}", "{assembly}.confident_cycs.fasta")),
               caption="../report/scapp.rst", category = "Plasmids",
               subcategory="SCAPP")
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

rule recycler:
    input:
        graph = opj("results", "assembly", "{assembly}", "assembly_graph.fastg.gz"),
        bam = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam"),
        bai = opj("results", "assembly", "{assembly}", "reads_pe_primary.sort.bam.bai"),
        kmer = opj("results", "assembly", "{assembly}", "kmer")
    output:
        touch(report(opj("results", "recycler", "{assembly}", "{assembly}.cycs.fasta"),
                     caption="../report/recycler.rst", category = "Plasmids",
                     subcategory="Recycler"))
    log:
        opj("results", "logs", "plasmids", "{assembly}.recycler.log")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    params:
        account=config["project"],
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