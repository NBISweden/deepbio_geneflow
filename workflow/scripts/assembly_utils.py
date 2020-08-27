#!/usr/bin/env python

import pandas as pd
from Bio.SeqIO import parse
from pathlib import Path

def metaspades_input(sm):
    """
    Generates fastq files to use as input for metaspades assembler

    :param sm: snakemake object
    :return:
    """
    from common import rename_records
    files = {"R1": [], "R2": []}
    assembly_dict = sm.params.assembly
    # Collect all files belonging to the assembly group
    for sample in assembly_dict.keys():
        for unit in assembly_dict[sample]:
            for pair in assembly_dict[sample][unit].keys():
                files[pair].append(
                    assembly_dict[sample][unit][pair][0])
    # Rename and concatenate reads (required for Metaspades)
    with open(sm.output.R1, 'w') as fh1, open(sm.output.R2, 'w') as fh2:
        for i, f in enumerate(files["R1"]):
            f2 = files["R2"][i]
            fh1 = rename_records(f, fh1, i)
            fh2 = rename_records(f2, fh2, i)


def main(sm):
    toolbox = {"generate_metaspades_input": metaspades_input}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)