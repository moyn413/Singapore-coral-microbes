#----------------------------------------------------------------
# nifH phyloseq processing file and script to make plots with nifH data. Starting file is a phyloseq file made using the dada2 pipeline.
# Treemap, deseq, and normalization were based on scripts from https://github.com/slimelab/Singapore-metabarcodes
#----------------------------------------------------------------

library(plyr) 
library(dplyr)
library(DESeq2)
library(ggplot2)
library(phyloseq)
library(decontam)
library(stringr)
library(seqinr)
library(DECIPHER)
library(dada2)
library(gridExtra)
library(tidyverse)
library(knitr)
library(treemapify)
library(RColorBrewer)
library(metagMisc)


#--------------------------------
# Load phyloseq file
#--------------------------------

ps <- readRDS("/Users/molly/Documents/GitHub/Singapore-coral-microbes/nifH/SGcoral_nifH_phyloseq_cluster.rds")   #read RDS phyloseq file
# dataset_path <- "~/GitHub/Singapore-coral-microbes/nifH"
# path_dataset <- function(file_name) str_c(dataset_path, file_name)
# treefasta_dir <-    path_dataset("/processed_nifH_ASVs/")


#--------------------------------
# Remove contaminants    
#--------------------------------

#Using decontam function, identify contaminants using the prevalence method with a threshold of 0.1. Note that additional SW sequences were added to "balance" the number of coral samples and avoid a bias in the decontam function. These samples are identified in the phyloseq as Type2 = decontam. 
#Negative controls are identified as "Control" in the phyloseq. 
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant)    #no contaminants identified

ps.noncontam <- prune_taxa(!contamdf.prev01$contaminant, ps) 

#Remove samples uses as negative controls and SW samples used to balance decontam
ps.noncontam_noC <-  ps.noncontam %>% subset_samples(Sample_or_Control != "Control") %>% subset_samples(Type2 != "decontam") 

#Rename main phyloseq
ps.all <- ps.noncontam_noC
ps.coral.all <-  ps.all %>% subset_samples(Type != "Seawater")
ps.coral.dna <-  ps.coral.all %>% subset_samples(NucleicType == "DNA")  
ps.coral.rna <-  ps.coral.all  %>% subset_samples(NucleicType == "RNA")
ps.SW <-  ps.all %>% subset_samples(Type == "Seawater")

# Remove samples with fewer than 20 reads 
ps.all.LR <- prune_samples(sample_sums(ps.all)>=20, ps.all)
ps.coral.all.LR <- prune_samples(sample_sums(ps.coral.all)>=20, ps.coral.all)
ps.coral.dna.LR <- prune_samples(sample_sums(ps.coral.dna)>=20, ps.coral.dna)
ps.coral.rna.LR <- prune_samples(sample_sums(ps.coral.rna)>=20, ps.coral.rna)
ps.SW.LR <- prune_samples(sample_sums(ps.SW)>=20, ps.SW)

# Remove OTUs occurring fewer than 10x
ps.all.LO = prune_taxa(taxa_sums(ps.all.LR)>=10, ps.all.LR)
ps.coral.all.LO <- prune_taxa(taxa_sums(ps.coral.all.LR)>=10, ps.coral.all.LR)
ps.coral.dna.LO <- prune_taxa(taxa_sums(ps.coral.dna.LR)>=10, ps.coral.dna.LR)
ps.coral.rna.LO <- prune_taxa(taxa_sums(ps.coral.rna.LR)>=10, ps.coral.rna.LR)
ps.SW.LO <- prune_taxa(taxa_sums(ps.SW.LR)>=10, ps.SW.LR) 

#------------------------------------------------------------------------------------------------
#Export sequences as fasta files and inspect in Geneious Prime
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.all.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "ALL_nifH_ASV_LO.fasta"), compress = FALSE, width = 20000)
# 
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.coral.dna.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "CORALnifH_ASV_DNA_LO.fasta"), compress = FALSE, width = 20000)
# 
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.coral.rna.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "CORALnifH_ASV_RNA_LO.fasta"), compress = FALSE, width = 20000)
# 
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.SW.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "CORALnifH_ASV_SW_LO.fasta"), compress = FALSE, width = 20000)

#Add Newick file made in Geneious Prime to phyloseq
all.dna.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/nifH/nifH_August_coral_SW_plus/dada2_phylum/trees/ALL_nifH_ASV_LO alignment FastTree Tree.newick")
coral.dna.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/nifH/nifH_August_coral_SW_plus/dada2_phylum/trees/CORALnifH_ASV_DNA_LO alignment FastTree Tree.newick")
coral.rna.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/nifH/nifH_August_coral_SW_plus/dada2_phylum/trees/CORALnifH_ASV_RNA_LO alignment FastTree Tree.newick")
SW.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/nifH/nifH_August_coral_SW_plus/dada2_phylum/trees/CORALnifH_ASV_SW_LO alignment FastTree Tree.newick")

#------------------------------------------------------------------------------------------------
ps.all.LO <- merge_phyloseq(ps.all.LO, all.dna.tree) 
ps.coral.dna.LO <- merge_phyloseq(ps.coral.dna.LO, coral.dna.tree) 
ps.coral.rna.LO <- merge_phyloseq(ps.coral.rna.LO, coral.rna.tree) 
ps.SW.LO <- merge_phyloseq(ps.SW.LO, SW.tree) 


# Remove OTUs identified in Geneious Prime as being erroneous
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

# remove outlier OTUs
badTaxa = c("otu2803", "otu1600", "otu3171")
ps.all.LO = pop_taxa(ps.all.LO, badTaxa)
ps.coral.dna.LO = pop_taxa(ps.coral.dna.LO, badTaxa)
ps.coral.rna.LO = pop_taxa(ps.coral.rna.LO, badTaxa)
ps.SW.LO = pop_taxa(ps.SW.LO, badTaxa)


#------------------------------------------------------------------------------------------------
# Make nifH phyloseqs for clusters I - III only & normalize
#------------------------------------------------------------------------------------------------

#Make phyloseq with ASVs from Cluster I - III only (clusters capable of fixing nitrogen)
ps.nfix <- subset_taxa(ps.all.LO, Cluster !="Cluster V chL/bchL" & Cluster !="Cluster V chlX/bchX" & Cluster != "Cluster VI"  & Cluster != "Cluster IV" & Cluster != "unknown")

#Make nfix phyloseq for coral DNA & RNA, coral DNA only, coral RNA only, and seawater 
ps.nfix.coral <-ps.nfix %>% subset_samples(Type != "Seawater")
ps.nfix.coral.dna <- ps.nfix.coral%>% subset_samples(NucleicType == "DNA")  
ps.nfix.coral.rna <- ps.nfix.coral%>% subset_samples(NucleicType == "RNA")  
ps.nfix.SW <-  ps.nfix %>% subset_samples(Type == "Seawater")

#Reprune to remove any samples with 0s
ps.nfix <- prune_samples(sample_sums(ps.nfix)>=1, ps.nfix)
ps.nfix.coral <- prune_samples(sample_sums(ps.nfix.coral)>=1, ps.nfix.coral)
ps.nfix.coral.dna <- prune_samples(sample_sums(ps.nfix.coral.dna)>=1, ps.nfix.coral.dna)
ps.nfix.coral.rna <- prune_samples(sample_sums(ps.nfix.coral.rna)>=1, ps.nfix.coral.rna)
ps.nfix.SW <- prune_samples(sample_sums(ps.nfix.SW)>=1, ps.nfix.SW)

#Some otus will be 0 after this selection, so re-prune to remove zeros 
ps.nfix <- prune_taxa(taxa_sums(ps.nfix)>=1, ps.nfix)
ps.nfix.coral <- prune_taxa(taxa_sums(ps.nfix.coral)>=1, ps.nfix.coral)
ps.nfix.coral.dna <- prune_taxa(taxa_sums(ps.nfix.coral.dna)>=1, ps.nfix.coral.dna)
ps.nfix.coral.rna <- prune_taxa(taxa_sums(ps.nfix.coral.rna)>=1, ps.nfix.coral.rna)
ps.nfix.SW <- prune_taxa(taxa_sums(ps.nfix.SW)>=1, ps.nfix.SW)



#-----------Normalize functions------
ps_normalize_median <- function (ps, title) {
  ps_median = median(sample_sums(ps))
  cat(sprintf("\nThe median number of reads used for normalization of %s is  %.0f", title, ps_median))
  normalize_median = function(x, t=ps_median) (if(sum(x) > 0){ t * (x / sum(x))} else {x})
  ps = transform_sample_counts(ps, normalize_median)
  cat(str_c("\nPhyloseq ",title, "\n========== \n") )
  print(ps)
}

ps_abundant <- function(ps,contrib_min=0.01, title){
  total_per_sample <- max(sample_sums(ps))
  ps <- filter_taxa(ps, function(x) sum(x > total_per_sample*contrib_min) > 0, TRUE)
  ps <- ps_normalize_median(ps, title)
}


#-----------Normalize phyloseq with all clusters ---------------------------------

# Normalize phyloseq files with low abundant OTUs (for diversity plots)
ps.all.LR.norm = ps_normalize_median(ps.all.LR, "nifH_all_LR_norm")
ps.coral.all.LR.norm = ps_normalize_median(ps.coral.all.LR, "Coral_nifH_all_LR_norm")
ps.coral.dna.LR.norm = ps_normalize_median(ps.coral.dna.LR, "Coral_nifH_dna_LR_norm")
ps.coral.rna.LR.norm = ps_normalize_median(ps.coral.rna.LR, "Coral_nifH_rna_LR_norm")
ps.SW.LR.norm = ps_normalize_median(ps.SW.LR, "SW_nifH_LR_norm")

# Normalize phyloseq file without low abundant OTUs
ps.all.norm = ps_normalize_median(ps.all.LO, "nifH_all_norm")
ps.coral.all.norm = ps_normalize_median(ps.coral.all.LO, "Coral_nifH_all_norm")
ps.coral.dna.norm = ps_normalize_median(ps.coral.dna.LO, "Coral_nifH_dna_norm")
ps.coral.rna.norm = ps_normalize_median(ps.coral.rna.LO, "Coral_nifH_rna_norm")
ps.SW.norm = ps_normalize_median(ps.SW.LO, "SW_nifH_norm")

# Remove non-abundant taxa and normalize again
ps.all.abund = ps_abundant(ps.all.norm,contrib_min=0.01, "nifH_all_abund")  
ps.coral.all.abund = ps_abundant(ps.coral.all.norm,contrib_min=0.01, "Coral_nifH_all_abund")
ps.coral.dna.abund = ps_abundant(ps.coral.dna.norm,contrib_min=0.01, "Coral_nifH_dna_abund")
ps.coral.rna.abund = ps_abundant(ps.coral.rna.norm,contrib_min=0.01, "Coral_nifH_rna_abund")
ps.SW.abund = ps_abundant(ps.SW.norm,contrib_min=0.01, "SW_nifH_abund") 


#-----------Normalize phyloseq that have clusters I to III only ---------------------------------
ps.nfix.norm  = ps_normalize_median(ps.nfix, "nifH_clusters_norm")
ps.nfix.coral.norm  = ps_normalize_median(ps.nfix.coral, "Coral_nifH_clusters_norm")
ps.nfix.coral.dna.norm = ps_normalize_median(ps.nfix.coral.dna, "Coral_nifH_clusters_dna_norm")
ps.nfix.coral.rna.norm = ps_normalize_median(ps.nfix.coral.rna, "Coral_nifH_clusters_rna_norm")
ps.nfix.SW.norm = ps_normalize_median(ps.nfix.SW, "SW_nifH_clusters_norm")





##---------------------------------------------------------------------
##----------MAIN MANUSCRIPT FIGURES-------------------------------
##---------------------------------------------------------------------


#------------------------------------------------------------------------------------------------
# Figure 5. Treemap of Cluster I--III ASVs of the \textit{nifH} DNA-based community for each site
#------------------------------------------------------------------------------------------------


#-----Functions for treemap--------#
ps_to_long <- function(ps) {
  otu_df <- data.frame(otu_table(ps)) %>% 
    rownames_to_column(var = "otu_id")
  taxo_df <- data.frame(tax_table(ps))%>% 
    rownames_to_column(var = "otu_id")
  otu_df <- left_join(taxo_df, otu_df)
  otu_df <- gather(otu_df, "sample","n_seq" ,contains("X")) # All samples contain X
  metadata_df <- data.frame(sample_data(ps)) %>% 
    rownames_to_column(var = "sample")
  otu_df <- left_join(otu_df, metadata_df)
}



# Function for Treemap with custom color palette 
treemap_gg_dv2 <- function(df, group1, group2, title, c.palette) {
  group1 <- enquo(group1)
  group2 <- enquo(group2)
  
  df <- df %>% group_by(!!group1, !!group2) %>%  summarise(n_seq=sum(n_seq))
  g_treemap <- ggplot(df, aes(area = n_seq, fill = !!group2, 
                              label = !!group2, subgroup = !!group1)) +
    ggtitle(title) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_subgroup_border() +
    treemapify::geom_treemap_text(colour = "black", place = "topleft", reflow = T, 
                                  padding.x =  grid::unit(3, "mm"), 
                                  padding.y = grid::unit(3, "mm") ) +
    treemapify::geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour =
                                             "white", fontface = "italic", min.size = 0) +
    scale_fill_manual(values=c.palette) +
    theme(legend.position="none", plot.title = element_text(size = 16, face = "bold"))
  print(g_treemap)
  return(g_treemap)
}



#-----Make ps_to_long for each species/type (tissue/skeleton/site) of Cluster I - III only--------#
#Platygyra Kusu tissue, replace NA with "unknown Order"
platyK.dna.tissue <- subset_samples(ps.nfix.coral.dna.norm, Species=="K_Platy" & Type == "Tissue")
platyK.dna2.tissue <- ps_to_long(platyK.dna.tissue)%>% replace_na(list(Order= "unknown Order"))

#Platygyra Hantu tissue, replace NA with "unknown Order"
platyH.dna.tissue <- subset_samples(ps.nfix.coral.dna.norm, Species=="H_Platy" & Type == "Tissue")
platyH.dna2.tissue <- ps_to_long(platyH.dna.tissue)%>% replace_na(list(Order= "unknown Order"))

#Goniopora Hantu tissue, replace NA with "unknown Order"
goni.dna.tissue <- subset_samples(ps.nfix.coral.dna.norm, Species=="H_Goni" & Type == "Tissue")
goni.dna2.tissue <- ps_to_long(goni.dna.tissue)%>% replace_na(list(Order= "unknown Order"))

#Pocillopora Kusu skeleton, replace NA with "unknown Order"
pocil.dna.skeleton <- subset_samples(ps.nfix.coral.dna.norm, Species=="K_Pocil" & Type == "Skeleton")
pocil.dna2.skeleton <- ps_to_long(pocil.dna.skeleton)%>% replace_na(list(Order= "unknown Order"))

#Platygyra Kusu skeleton, replace NA with "unknown Order"
platyK.dna.skeleton <- subset_samples(ps.nfix.coral.dna.norm, Species=="K_Platy" & Type == "Skeleton")
platyK.dna2.skeleton <- ps_to_long(platyK.dna.skeleton)%>% replace_na(list(Order= "unknown Order"))

#Platygyra Hantu skeleton, replace NA with "unknown Order"
platyH.dna.skeleton <- subset_samples(ps.nfix.coral.dna.norm, Species=="H_Platy" & Type == "Skeleton")
platyH.dna2.skeleton <- ps_to_long(platyH.dna.skeleton)%>% replace_na(list(Order= "unknown Order"))

#Goniopora Hantu skeleton, replace NA with "unknown Order"
goni.dna.skeleton <- subset_samples(ps.nfix.coral.dna.norm, Species=="H_Goni" & Type == "Skeleton")
goni.dna2.skeleton <- ps_to_long(goni.dna.skeleton)%>% replace_na(list(Order= "unknown Order"))



#-----Make color palette for treemap--------#
getPalette = colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))) 
orderList = unique(tax_table(ps.nfix.coral.norm)[,"Order"])
orderPalette = getPalette(length(orderList))
names(orderPalette) = orderList


#-----Make treemap for each coral and type (tissue/skeleton)--------#

# (Pocillopora tissue has no cluster I - III ASVs) so it is not processed here
platyK.tree.dna.tissue <- treemap_gg_dv2(platyK.dna2.tissue, Cluster, Order, "Platygyra Kusu - Tissue DNA samples", orderPalette)
platyH.tree.dna.tissue <- treemap_gg_dv2(platyH.dna2.tissue, Cluster, Order, "Platygyra Hantu - Tissue DNA samples", orderPalette)
goni.tree.dna.tissue <- treemap_gg_dv2(goni.dna2.tissue, Cluster, Order, "Goniopora Hantu - Tissue DNA samples", orderPalette)

pocil.tree.dna.skeleton <- treemap_gg_dv2(pocil.dna2.skeleton, Cluster, Order, "Pocillopora Kusu - Skeleton DNA samples", orderPalette)
platyK.tree.dna.skeleton <- treemap_gg_dv2(platyK.dna2.skeleton, Cluster, Order, "Platygyra Kusu - Skeleton DNA samples", orderPalette)
platyH.tree.dna.skeleton <- treemap_gg_dv2(platyH.dna2.skeleton, Cluster, Order, "Platygyra Hantu - Skeleton DNA samples", orderPalette)
goni.tree.dna.skeleton <- treemap_gg_dv2(goni.dna2.skeleton, Cluster, Order, "Goniopora Hantu - Skeleton DNA samples", orderPalette)

#Arrange plots in grid
treemap.grid <- grid.arrange(goni.tree.dna.tissue, platyH.tree.dna.tissue, platyK.tree.dna.tissue, pocil.tree.dna.skeleton,
                             goni.tree.dna.skeleton, platyH.tree.dna.skeleton, platyK.tree.dna.skeleton, pocil.tree.dna.skeleton, nrow = 2)

# ggsave("nifH_I-III_DNA_species_tissue_skeleton_order_cluster_updated.pdf", plot = treemap.grid,
#        path = "/path",
#        width = 30,
#        height = 18)


#------------------------------------------------------------------------------------------------
# Figure 6. nifH deseq between the skeleton and tissue (Clusters I - III only)
#------------------------------------------------------------------------------------------------


ps_do_deseq <- function(deseq_data, a){
  
  diagdds = phyloseq_to_deseq2(deseq_data, ~Type2)   # Choose skeleton or tissue or DNA/RNA
  
  #---------Estimate geometric mean with zeros
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq(diagdds, fitType="local")
  
  #---------Finish loading deseq file
  
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = a
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(deseq_data)[rownames(sigtab), ], "matrix"))  ##change ps name
  head(sigtab)
  dim(sigtab)
  theme_bw()
  #---------
  
  
  #--------- Select color palette
  
  getPalette = c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"), brewer.pal(9, "Set1"), brewer.pal(8, "Spectral"), brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
  
  
  #--------- Plot
  # Class & Order plots 
  x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
  g <- ggplot(sigtab, aes(x=Order, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+   
    ylab("log2 fold change - Skeleton <--> Tissue") +
    scale_colour_manual(values=getPalette) +  coord_flip() + theme_bw() +  
    theme(axis.text.x = element_text(size = 10),
          axis.title.x =element_text(size=10),
          axis.text.y = element_text(size = 10 ),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10)) + 
          scale_y_continuous(limits = c(-7,7))   # Changed scale here to balance plot y-axis around 0, but one outlier does not fit in this new scale. Saved outlier information and inserted seperately using Illustrator for final plot. To see all values, remove the scale_y_continuous term here.
  return(g)
}




#Make deseq comparing tissue and skeleton nifH for Clusters I - III only. Because Pocillopora tissue has no Cluster I - III nifH in the tissue, 
# exlude Pocillopora from the tissue/skeleton comparison

deseq.subset <- subset_samples(ps.nfix.coral.dna, Species != "K_Pocil")
nfix.deseq <- ps_do_deseq(deseq.subset, 0.8) + ggtitle("all samples clusters I - III DNA   alpha= 0.8")


# ggsave("deseq_all_nifH_I-III_DNA_order_nopocil.pdf", plot = nfix.deseq,
#        path = "/path",
#        width = 10,
#        height = 7)


##---------------------------------------------------------------------
##----------SUPPLEMENTAL FIGURES-------------------------------
##---------------------------------------------------------------------

##---------------------------------------------------------------------
# Figure S14. nifH Cluster abundance (all clusters)
##---------------------------------------------------------------------


#Bar plot with facet wrap by clusters
dna.all.norm <- subset_samples(ps.all.LR.norm, NucleicType == "DNA")
top <- names(sort(taxa_sums(ps.all.LR.norm), decreasing=TRUE))
ps.top <- transform_sample_counts(dna.all.norm, function(OTU) OTU/sum(OTU))
sup.bar <- plot_bar(ps.top, fill="Cluster", x="Bars3") + facet_wrap(~Cluster, scales="free_x")+
  theme(legend.key.size = unit(2,"line"))  + theme(legend.position="bottom") + theme_bw()+
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +  ggtitle("nifH DNA all") + 
  scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")) + scale_color_manual (values =c(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")))


# ggsave("nifH_cluster_barplots_all.pdf", plot =sup.bar, path = "/path",
#        width = 15,
#        height = 10)


##---------------------------------------------------------------------
# Figure S15.  nifH RNA-based Clusters I--III
##---------------------------------------------------------------------
getPalette = colorRampPalette(c(brewer.pal(8,"Dark2"),brewer.pal(9, "Set1"),brewer.pal(12, "Set3"),brewer.pal(8,"Spectral")))
classList = unique(tax_table(ps.nfix.coral.norm)[,"Class"])
classPalette = getPalette(length(classList))
names(classPalette) = classList


top <- names(sort(taxa_sums(ps.nfix.coral.norm), decreasing=TRUE))
ps.top <- transform_sample_counts(ps.nfix.coral.norm, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(top, ps.top)
sup.bar2 <- plot_bar(ps.top, fill="Class", x="Sample") + facet_wrap(NucleicType~Cluster, scales="free_x")+
  theme(legend.key.size = unit(0.9,"line"))  + theme(legend.position="bottom") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  theme(strip.text.x = element_text(size = 12))+ scale_color_manual(values=classPalette) + scale_fill_manual(values= classPalette)


# ggsave("nifH_clusterI-III_barplots_Class.pdf", plot =sup.bar2, path = "/path",
#        width = 10,
#        height = 8)


##---------------------------------------------------------------------
# Figure S17. nifH RNA:DNA ratios of Clusters I - III only
##---------------------------------------------------------------------

#Select samples for which both RNA and DNA data are available
ratio.nfix.coral <- ps.nfix.coral.norm%>% subset_samples(Name =="14S"| Name =="15S"| Name =="15T"| Name == "2S"| Name == "6T"| Name == "7T"| Name == "8T")

# Make custom color palette for the order level
getPalette = colorRampPalette(c(brewer.pal(12, "Paired"), brewer.pal(8, "Set1"), brewer.pal(8, "Pastel2"), brewer.pal(8, "Accent"))) 
orderList = unique(tax_table(ratio.nfix.coral)[,"Order"])
orderPalette = getPalette(length(orderList))
names(orderPalette) = orderList

#Merge replicate samples
merged.nfix <- merge_samples(ratio.nfix.coral, "NucleicType")

#Normalize
nfix.norm = ps_normalize_median(merged.nfix, "merged_coral_norm")

#Convert phyloseq to dataframe
nfix.ASV.df <-phyloseq_to_df(nfix.norm, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance")

#Select only ASVs present in both RNA and DNA samples
nfix.ASV.df.select <-nfix.ASV.df %>% filter(RNA > 0 & DNA >0)

s17 <- ggplot(nfix.ASV.df.select, aes(x= DNA, y= RNA, colour=Class, label=Order))+
  geom_point(size = 3) + geom_text(aes(label=Order),hjust=0, vjust=0) +  geom_abline(intercept = 0) +
  theme_bw() +  ggtitle("RNA:DNA nifH")+
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10')+
  annotation_logticks(sides="lb")+  expand_limits(x = 1000 ,y = 1000)  +  scale_color_manual(values=orderPalette)

# ggsave("nifH_RNA_DNA_byOTU.pdf", plot = s17, path = "/path",
#        width = 8,
#        height = 6)
