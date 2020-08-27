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
    import shutil
    files = []
    assembly_dict = sm.params.assembly
    # Collect all files belonging to the assembly group
    for sample in sorted(assembly_dict.keys()):
        for unit in sorted(assembly_dict[sample]):
            f = assembly_dict[sample][unit][sm.params.R][0]
            files.append(f)
    # Rename and concatenate reads (required for Metaspades)
    with open(sm.params.tmp_out, 'w') as fh:
        for i, f in enumerate(files):
            fh = rename_records(f, fh, i)
    shutil.move(sm.params.tmp_out, sm.output.fq)


def main(sm):
    toolbox = {"generate_metaspades_input": metaspades_input}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)