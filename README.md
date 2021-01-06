![CI](https://github.com/NBISweden/deepbio_geneflow/workflows/CI/badge.svg)

# Deep Biosphere Geneflow project
Analyzing gene flow in the deep biosphere via "multi-omics" integrated analysis.

## Repository organization
`config/` : configuration files

`data/` : raw data goes here

`resources/` : databases, references etc.

`workflow/` : main directory for the snakemake workflow

## Data locations

`data/testdata/` : 

This folder contains two small test datasets:

- **test**

This dataset is comprised of 25k paired-end reads 
from the [SCAPP](https://github.com/Shamir-Lab/SCAPP) test data + 25k paired-end 
reads each from three plasmids:
- CEX4 plasmid pCEX4 (LC556220.1, *Enterobacter cloacae*)
- unnamed plasmid (NC_012780.1, *[Eubacterium] eligens ATCC 27750*)
- plasmid pBPSE01 (NZ_KF418775.1, *Burkholderia pseudomallei strain MSHR1950*)
generated using `randomreads.sh` from `bbmap`.

- **mock**

comprised of 50k reads from the `minced` testdata of the *Aquifex aeolicus* VF5 
genome (generated with `randomreads.sh`) + 50k reads subsampled from a 
[synthetic mock]() metagenome (using `seqtk`).
 
 ## Tools and outline
 ### Plasmids

- [SCAPP](https://github.com/Shamir-Lab/SCAPP) ([Pellow et al 2020](https://www.biorxiv.org/content/10.1101/2020.01.12.903252v3))
  
  SCAPP (Sequence Contents-Aware Plasmid Peeler) also uses the assembly graph output together with
  plasmid-specific gene detection and scoring of putative plasmids.

- [Recycler](https://github.com/Shamir-Lab/Recycler) ([Rozov et al 2016](https://academic.oup.com/bioinformatics/article/33/4/475/2623362))
  
  Recycler also uses assembly graphs + read-mapping and coverage information to identify 
  putative plasmids but does not use plasmid-specific genes.
  
- [PPR-Meta](https://github.com/zhenchengfang/PPR-Meta) ([Fang et al 2019](https://academic.oup.com/gigascience/article/8/6/giz066/5521157))
  
  PPR-Meta uses deep learning methods to identify phages and plasmid sequences in 
  metagenomic assemblies.

### Phages
- [Virsorter](https://github.com/simroux/VirSorter) ([Roux et al 2015](https://peerj.com/articles/985/))

  Virsorter uses curated protein databases of viral genes which are queried with
  `hmmsearch` or `blastp` using genes predicted on contigs as queries. Metrics
  are then computed using sliding windows to look for enrichment of 'viral' 
  signals, followed by calculation of a significance score. It also includes a 
  step to identify circular sequences as an initial step. The author notes that:
  > for fragmented genomes, *category 3* predictions help recover more viral 
  >sequences, but do so at the cost of increased false-positives.

- [MARVEL](https://github.com/LaboratorioBioinformatica/MARVEL) ([Amgarten et al 2018](https://www.frontiersin.org/articles/10.3389/fgene.2018.00304/full)) 
  
  Marvel uses a random forest classifier to classify metagenomic bins as 
  phage/bacteria based on features such as 1) gene density, 2) strand shifts
  and 3) fraction of hits to [pVOG database](http://dmk-brain.ecn.uiowa.edu/pVOGs/).
  The classifier was trained on 1,247 phage and 1,029 bacterial genomes. The 
  authors note that:
  >MARVEL has high F1 scores and accuracy for all bin lengths analyzed, but 
  >especially for bins composed of contigs 4 kbp long and longer.
 
- [VIBRANT](https://github.com/AnantharamanLab/VIBRANT/) ([Kieft et al 2020](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0))

  VIBRANT uses machine learning (neural networks) based on annotation
  metrics derived from HMM searches against KEGG, PFAM AND [VOG](http://vogdb.org/).  

### CRISPRs
#### [Minced](https://github.com/ctSkennerton/minced)

### Genome Islands
- Benchmarking shows that assembly + binning underperforms for plasmids and GIs: 
[Maguire et al 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.997171v2.abstract)