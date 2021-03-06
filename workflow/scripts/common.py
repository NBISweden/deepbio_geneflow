#!/usr/bin/env python

from os.path import join as opj


def parse_samples(df):
    assemblies = {}
    samples = {}
    df.fillna("", inplace=True)

    for i in list(df.index):
        # Add sample to dict
        sample = df.iloc[i]["sample"]
        if sample not in samples.keys():
            samples[sample] = {}
        # Add unit to dict
        unit = str(df.iloc[i]["unit"])
        if unit not in samples[sample].keys():
            samples[sample][unit] = {}
        samples[sample][unit]["R1"] = df.iloc[i]["fq1"]
        samples[sample][unit]["R2"] = df.iloc[i]["fq2"]

        # Initiate keys for all assembly group values
        assem_list = df.iloc[i]["assembly"].split(",")
        assem_list = [a for a in assem_list if a != ""]
        for a in assem_list:
            if a not in assemblies.keys():
                assemblies[a] = {}
            if sample not in assemblies[a].keys():
                assemblies[a][sample] = {unit: {}}
            if unit not in assemblies[a][sample].keys():
                assemblies[a][sample][unit] = {}
            assemblies[a][sample][unit]["R1"] = [samples[sample][unit]["R1"]]
            assemblies[a][sample][unit]["R2"] = [samples[sample][unit]["R2"]]
    return samples, assemblies


def get_assembly_files(assembly_dict, R):
    files = []
    for sample in sorted(assembly_dict.keys()):
        for unit in sorted(assembly_dict[sample].keys()):
            files.append(assembly_dict[sample][unit][R][0])
    return files


def rename_records(f, fh, i):
    """
    Prepends a number to read ids

    :param f: Input fastq file (gzipped)
    :param fh: Output filehandle
    :param i: File index to prepend to reads
    :return: Output filehandle
    """
    from Bio import SeqIO
    import gzip as gz
    for record in SeqIO.parse(gz.open(f, 'rt'), 'fastq'):
        record.id = "{}_{}".format(i, record.id)
        SeqIO.write(record, fh, "fastq")
    return fh
