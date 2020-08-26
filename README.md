# Deep Biosphere Geneflow project
Analyzing gene flow in the deep biosphere via "multi-omics" integrated analysis.

## Repository organization
`config/` : configuration files

`data/` : raw data goes here

`resources/` : databases, references etc.

`workflow/` : main directory for the snakemake workflow

## Data locations

 ## Tools and outline
 ### Plasmids

- [MetaplasmidSPADES](https://github.com/ablab/spades/tree/metaplasmid_3.13.0) ([Antipov et al 2019](https://genome.cshlp.org/content/29/6/961))
- [SCAPP](https://github.com/Shamir-Lab/SCAPP) ([Pellow et al 2020](https://www.biorxiv.org/content/10.1101/2020.01.12.903252v3))
- [Recycler](https://github.com/Shamir-Lab/Recycler) ([Rozov et al 2016](https://academic.oup.com/bioinformatics/article/33/4/475/2623362))

### Phages
- [Virsorter](https://github.com/simroux/VirSorter) ([Roux et al 2015](https://peerj.com/articles/985/))
- [MetaviralSPADES](https://github.com/ablab/spades/tree/metaviral_publication) ([Antipov et al 2020](https://academic.oup.com/bioinformatics/article-abstract/36/14/4126/5837667)) 

### CRISPR
- [Minced](https://github.com/ctSkennerton/minced)

### Genome Islands
- Benchmarking shows that assembly + binning underperforms for plasmids and GIs: 
[Maguire et al 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.997171v2.abstract)