# Singapore-coral-microbes

Script and data for: 

Moynihan, M.A., Goodkin, N.F., Morgan, K.M., Kho, P.Y.Y, Lopes dos Santos, A., Lauro, F.M., Baker, D.M, and Martin, P. Coral-associated nitrogen fixation rates and diazotrophic diversity on a nutrient-replete equatorial reef. (in prep)

## Directories

### nifH

https://github.com/moyn413/Singapore-coral-microbes/tree/main/nifH/

* Newick and alignment files used to make nifH phylogeny based on protein alignment of select nifH sequences and reference sequences. Reference sequences in the tree are based on Meheust et al. 2020. ([SG_nifH_protein_tree_FINAL.newick](https://github.com/moyn413/Singapore-coral-microbes/tree/main/nifH/TREE-nifH_protein_phylogeny/SG_nifH_protein_tree_FINAL.newick ) & [SG_nifH_protein_alignment.fasta](https://github.com/moyn413/Singapore-coral-microbes/tree/main/nifH/TREE-nifH_protein_phylogeny/SG_nifH_protein_alignment.fasta))
* fasta file with all nifH sequences ([SG_nifH_ASV.fasta](https://github.com/moyn413/Singapore-coral-microbes/tree/main/nifH/blast/SG_nifH_ASV.fasta))
* ASV table, where nifH clusters, where manually assigned ([nifH_dada2_ASV_table_wClusters.tsv](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/nifH_dada2_ASV_table_wClusters.tsv))
* Phyloseq used to make figures ([SG_nifH_phyloseq_cluster.rds](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/SG_nifH_phyloseq_cluster.rds))
* R script used to process phyloseq and make figures ([nifH_plots.R.rds](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/nifH_plots.R))
* Files generated & used during phyloseq processing ([processed_nifH_ASVs](https://github.com/moyn413/Singapore-coral-microbes/blob/main/nifH/processed_nifH_ASVs))

### 16S

* fasta file with all 16S sequences ([SGcoral_16S_ASV.fasta](https://github.com/moyn413/Singapore-coral-microbes/blob/master/16S/blast/SGcoral_16S_ASV.fasta))
* 16S ASV table ([SGcoral_16S_dada2_ASV_table.tsv](https://github.com/moyn413/Singapore-coral-microbes/blob/master/16S/SGcoral_16S_dada2_ASV_table.tsv))
* Phyloseq used to make figures ([SGcoral_16S_phyloseq.rds](https://github.com/moyn413/Singapore-coral-microbes/blob/master/16S/SGcoral_16S_phyloseq.rds))
Rmd
* R script used to process phyloseq and make figures ([SGcoral_16S_plots.R](https://github.com/moyn413/Singapore-coral-microbes/blob/master/16S/SGcoral_16S_plots.R))
* Files generated & used during phyloseq processing ([processed_16S_ASVs](https://github.com/moyn413/Singapore-coral-microbes/blob/master/16S/processed_16S_ASVs))
* Plastid 16S sequences, phyloseq, ASV table, and script ([Plastids](https://github.com/moyn413/Singapore-coral-microbes/blob/master/16S/Plastids))

### 18S 

* phyloseq.Rmd
* ASV_table.xlsx
* .R-script for making nifH plots in paper



## References
Méheust R, Castelle CJ, Carnevali PBM, Farag IF, He C, Chen LX, et al. Groundwater Elusimicrobia are metabolically diverse compared to gut microbiome Elusimicrobia and some have a novel nitrogenase paralog. The ISME Journal. 2020;p. 1–16.
