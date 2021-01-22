# Singapore-coral-microbes

Script and data for: 

Moynihan, M.A., Goodkin, N.F., Morgan, K.M., Kho, P.Y.Y, Lopes dos Santos, A., Lauro, F.M., Baker, D.M, and Martin, P. Coral-associated nitrogen fixation rates and diazotrophic diversity on a nutrient-replete equatorial reef. (in revision)

## Directories

### nifH

* Newick and alignment files used to make nifH phylogeny based on protein alignment of select nifH sequences and reference sequences. Reference sequences in the tree are based on Meheust et al. 2020. ([SGcoral_nifH_protein_tree.newick](https://github.com/moyn413/Singapore-coral-microbes/tree/main/nifH/nifH_protein_phylogeny/SGcoral_nifH_protein_tree.newick) & [SG_nifH_protein_alignment.fasta](https://github.com/moyn413/Singapore-coral-microbes/tree/main/nifH/nifH_protein_phylogeny/SGcoral_nifH_protein_alignment.fasta))
* fasta file with all nifH sequences ([SGcoral_nifH_ASV.fasta](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/SG_nifH_ASV.fasta))
* ASV table, where nifH clusters, where manually assigned ([SGcoral_nifH_dada2_ASV_table_wClusters.tsv](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/SGcoral_nifH_dada2_ASV_table_wClusters.tsv))
* Phyloseq used to make figures ([SGcoral_nifH_phyloseq_cluster.rds](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/SGcoral_nifH_phyloseq_cluster.rds))
* R script used to process phyloseq and make figures ([SGcoral_nifH_plots.R](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/SGcoral_nifH_plots.R))
* Files generated & used during phyloseq processing ([phyloseq tree files](https://github.com/moyn413/Singapore-coral-microbes/tree/main/nifH/phyloseq%20tree%20files))

### 16S

* fasta file with all 16S sequences ([SGcoral_16S_ASV.fasta](https://github.com/moyn413/Singapore-coral-microbes/blob/main/16S/SGcoral_16S_ASV.fasta))
* 16S ASV table ([SGcoral_16S_dada2_ASV_table.tsv](https://github.com/moyn413/Singapore-coral-microbes/blob/main/16S/SGcoral_16S_dada2_ASV_table.tsv))
* Phyloseq used to make figures ([SGcoral_16S_phyloseq.rds](https://github.com/moyn413/Singapore-coral-microbes/blob/main/16S/SGcoral_16S_phyloseq.rds))
* R script used to process phyloseq and make figures ([SGcoral_16S_plots.R](https://github.com/moyn413/Singapore-coral-microbes/blob/main/16S/SGcoral_16S_plots.R))
* Files generated & used during phyloseq processing ([phyloseq tree files](https://github.com/moyn413/Singapore-coral-microbes/tree/main/16S/phyloseq%20tree%20files))
* Plastid 16S sequences, phyloseq, ASV table, and script ([Plastids](https://github.com/moyn413/Singapore-coral-microbes/tree/main/16S/Plastids))

### 18S 

* fasta file with all 18S sequences ([SGcoral_18S_ASV.fasta](https://github.com/moyn413/Singapore-coral-microbes/blob/main/18S/SGcoral_18S_ASV.fasta))
* 18S ASV table ([SGcoral_18S_dada2_ASV_table.tsv](https://github.com/moyn413/Singapore-coral-microbes/blob/main/18S/SGcoral_18S_dada2_ASV_table.tsv))
* Phyloseq used to make figures ([SGcoral_18S_phyloseq.rds](https://github.com/moyn413/Singapore-coral-microbes/blob/main/18S/SGcoral_18S_phyloseq.rds))
* R script used to process phyloseq and make figures ([SGcoral_18S_plots.R](https://github.com/moyn413/Singapore-coral-microbes/blob/main/18S/SGcoral_18S_plots.R))



#### References
Méheust R, Castelle CJ, Carnevali PBM, Farag IF, He C, Chen LX, et al. Groundwater Elusimicrobia are metabolically diverse compared to gut microbiome Elusimicrobia and some have a novel nitrogenase paralog. The ISME Journal. 2020;p. 1–16.

Chénard, C., Wijaya, W., Vaulot, D. et al. Temporal and spatial dynamics of Bacteria, Archaea and protists in equatorial coastal waters. Sci Rep 9, 16390 (2019). 
