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

        # Set staged file paths
        r1 = opj("data", "stage", "{}_{}_R1.fastq.gz".format(sample, unit))
        r2 = opj("data", "stage", "{}_{}_R2.fastq.gz".format(sample, unit))

        # Initiate keys for all assembly group values
        assem_list = df.iloc[i]["assembly"].split(",")
        assem_list = [a for a in assem_list if a != ""]
        for a in assem_list:
            if a not in assemblies.keys():
                assemblies[a] = {}
            if sample not in assemblies[a].keys():
                assemblies[a][sample] = {unit: {}}
            assemblies[a][sample][unit]["R1"] = [r1]
            assemblies[a][sample][unit]["R2"] = [r2]
    return samples, assemblies


def get_assembly_files(assembly_dict):
    files = []
    for sample in assembly_dict.keys():
        for unit in assembly_dict[sample].keys():
            for pair in assembly_dict[sample][unit].keys():
                files.append(assembly_dict[sample][unit][pair][0])
    return files