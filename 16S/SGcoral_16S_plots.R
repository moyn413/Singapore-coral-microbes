#----------------------------------------------------------------
# 16S phyloseq processing file and script to make plots with 16S data. Starting file is a phyloseq file made using the dada2 pipeline.
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
library(vegan)

#--------------------------------
# Load phyloseq file
#--------------------------------

ps <- readRDS("/GitHub/Singapore-coral-microbes/16S/SGcoral_16S_phyloseq.rds")   #read RDS phyloseq file

dataset_path <- "/GitHub/Singapore-coral-microbes/16S/"
path_dataset <- function(file_name) str_c(dataset_path, file_name)
#dada2_dir <-  path_dataset("/16S/")
treefasta_dir <-path_dataset("/16S/phyloseq tree files/")

#--------------------------------
# Remove contaminants    
#--------------------------------

#Using decontam function, identify contaminants using the prevalence method with a threshold of 0.1. Note that additional SW sequences were added to "balance" the number of coral samples and avoid a bias in the decontam function. These samples are identified in the phyloseq as decontam. 
#Negative controls are identified as "Control" in the phyloseq. 
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant)    # 22 contaminants identified

# remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev01$contaminant, ps) 

#Remove samples uses as negative controls and SW samples used to balance decontam
ps.augustsamples <-  ps.noncontam %>% subset_samples(Sample_or_Control != "Control") %>% subset_samples(NucleicType != "decontam") 

#Remove eukaryote and chloroplast sequences. Subset taxa will remove NA so need to add exception for "or is NA"
ps.august_noeuk <- subset_taxa(ps.augustsamples, Domain != "Eukaryota"| is.na(Domain)) #ntaxa= 15720
ps.august_noeuk_noC <- subset_taxa(ps.august_noeuk, Order != "Chloroplast"| is.na(Order)) # ntaxa=15151
ps.august_noeuk_noC_noM <- subset_taxa(ps.august_noeuk_noC, Family != "Mitochondria"| is.na(Family)) #ntaxa = 14947

# Make separate phyloseqs for different sample types
ps.all <- ps.august_noeuk_noC_noM
ps.coral.all <-  ps.all %>% subset_samples(Type != "Seawater")
ps.coral.dna <-  ps.all %>% subset_samples(Type != "Seawater") %>% subset_samples(NucleicType == "DNA")
ps.coral.rna <-  ps.all %>% subset_samples(Type != "Seawater") %>% subset_samples(NucleicType == "RNA")
ps.SW <-  ps.all %>% subset_samples(Type == "Seawater")

# Remove samples with fewer than 20 reads (if needed)
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
ps.SW.LO <- prune_taxa(taxa_sums(ps.SW.LR)>10, ps.SW.LR)  

#------------------------------------------------------------------------------------------------
# Make fasta files and export to make Fast Tree in Geneious prime and import back into R to add to phyloseq
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.all.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "ALL_16S_ASV_LO.fasta"), compress = FALSE, width = 20000)
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.coral.dna.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "CORAL16S_ASV_DNA_LO.fasta"), compress = FALSE, width = 20000)
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.coral.rna.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "CORAL16S_ASV_RNA_LO.fasta"), compress = FALSE, width = 20000)
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.SW.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "CORAL16S_ASV_SW_LO.fasta"), compress = FALSE, width = 20000)
# seqs <- Biostrings::DNAStringSet(getSequences(refseq(ps.SW.LO)))
# Biostrings::writeXStringSet(seqs, str_c(treefasta_dir, "CORAL16S_ASV_SW_LO.fasta"), compress = FALSE, width = 20000)

#Add Newick files made in Geneious Prime to phyloseq
all.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/16S/16S_August_coral_SW_plus_wneg_new_silva/dada2/trees/phylo/ALL_16S_ASV_LO alignment FastTree Tree.newick")
coral.dna.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/16S/16S_August_coral_SW_plus_wneg_new_silva/dada2/trees/phylo/CORAL16S_ASV_DNA_LO alignment FastTree Tree.newick")
coral.rna.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/16S/16S_August_coral_SW_plus_wneg_new_silva/dada2/trees/phylo/CORAL16S_ASV_RNA_LO alignment FastTree Tree.newick")
SW.tree <- ape::read.tree("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/16S/16S_August_coral_SW_plus_wneg_new_silva/dada2/trees/phylo/CORAL16S_ASV_SW_LO alignment FastTree Tree.newick")
#------------------------------------------------------------------------------------------------


#Add tree into phyloseq
ps.all.LO.2 <- merge_phyloseq(ps.all.LO,all.tree)
ps.coral.dna.LO2 <- merge_phyloseq(ps.coral.dna.LO, coral.dna.tree) 
ps.coral.rna.LO2 <- merge_phyloseq(ps.coral.rna.LO, coral.rna.tree) 
ps.SW.LO2 <- merge_phyloseq(ps.SW.LO, SW.tree) 


# Remove OTUs identified in Geneious Prime as being erroneous  
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

badTaxa = c("otu6504", "otu5256","otu7142", "otu3459")

ps.all.LO = pop_taxa(ps.all.LO.2, badTaxa)
ps.coral.dna.LO = pop_taxa(ps.coral.dna.LO2, badTaxa)
ps.coral.rna.LO = pop_taxa(ps.coral.rna.LO2, badTaxa)
ps.SW.LO = pop_taxa(ps.SW.LO2, badTaxa)



#------------------------------------------------------------------------------------------------

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
#------------------------------------------------------------------------------------------------

#-------------------------------
#  Normalize 
#-------------------------------

# Normalize phyloseq files with low abundant OTUs (for diversity plots)
ps.all.LR.norm = ps_normalize_median(ps.all.LR, "16S_all_LR_norm")
ps.coral.all.LR.norm = ps_normalize_median(ps.coral.all.LR, "Coral_16S_all_LR_norm")
ps.coral.dna.LR.norm = ps_normalize_median(ps.coral.dna.LR, "Coral_16S_dna_LR_norm")
ps.coral.rna.LR.norm = ps_normalize_median(ps.coral.rna.LR, "Coral_16S_rna_LR_norm")
ps.SW.LR.norm = ps_normalize_median(ps.SW.LR, "SW_16S_norm")


# Normalize phyloseq file without low abundant OTUs
ps.all.norm = ps_normalize_median(ps.all.LO, "16S_all_norm")
ps.coral.dna.norm = ps_normalize_median(ps.coral.dna.LO, "Coral_16S_LO_dna_norm")
ps.coral.rna.norm = ps_normalize_median(ps.coral.rna.LO, "Coral_LO_16S_rna_norm")
ps.SW.norm = ps_normalize_median(ps.SW.LO, "SW_LO_16S_norm")

# Remove non-abundant taxa and normalize again
ps.all.abund = ps_abundant(ps.all.norm,contrib_min=0.01, "16S_all_abund")  
ps.coral.dna.abund = ps_abundant(ps.coral.dna.norm,contrib_min=0.01, "Coral_16S_dna_abund") 
ps.coral.rna.abund = ps_abundant(ps.coral.rna.norm,contrib_min=0.01, "Coral_16S_rna_abund")
ps.SW.abund = ps_abundant(ps.SW.norm,contrib_min=0.01, "SW_16S_abund") 

#Make coral allphyloseq that has tree (from ps.all.LO)
ps.coral.all.LO2 <-  ps.all.LO %>% subset_samples(Type != "Seawater")
ps.coral.all.norm = ps_normalize_median(ps.coral.all.LO2, "16S_coral_all_norm")
ps.coral.all.abund = ps_abundant(ps.coral.all.norm,contrib_min=0.01, "Coral_16S_all_abund")


##---------------------------------------------------------------------
##----------STATISTICS-------------------------------
##---------------------------------------------------------------------

#-------------------------------
#  Alpha Diversity 
#-------------------------------

# Perform stats on all samples, not trimmed. Use this for stats table. 
ps.diversity = transform_sample_counts(ps.all.LR, ceiling)  
rich = estimate_richness(ps.diversity)      # For Table S7


ps.diversity.dna <- ps.diversity %>% subset_samples(NucleicType == "DNA")
ps.diversity.rna <- ps.diversity %>% subset_samples(NucleicType == "RNA")

plot_richness(ps.diversity.dna, measures=c("Observed","Shannon", "Simpson"), x="SpeciesType", color="Type") + geom_boxplot()
plot_richness(ps.diversity.rna, measures=c("Observed","Shannon", "Simpson"), x="SpeciesType", color="Type") + geom_boxplot()


#-------------------------------
#  DNA Tissue vs Skeleton
#-------------------------------

pocil <- subset_samples(ps.all.LR, Species2=="K_Pocil" & NucleicType == "DNA")
platyK <- subset_samples(ps.all.LR, Species2=="K_Platy"& NucleicType == "DNA")
platyH <- subset_samples(ps.all.LR, Species2=="H_Platy"& NucleicType == "DNA")
goni <- subset_samples(ps.all.LR, Species2=="H_Goni"& NucleicType == "DNA")

# Kusu Pocillopora tissue/skeleton DNA, p=0.1
rich = estimate_richness(pocil, measures=c("Shannon", "Simpson", "Chao1"))
t.test(rich$Shannon~ sample_data(pocil)$Type)

# Kusu Platygyra tissue/skeleton DNA, p=0.01
rich = estimate_richness(platyK, measures=c("Shannon", "Simpson", "Chao1"))
t.test(rich$Shannon~ sample_data(platyK)$Type)

# Hantu Platygyra tissue/skeleton DNA, p=0.4
rich = estimate_richness(platyH, measures=c("Shannon", "Simpson", "Chao1"))
t.test(rich$Shannon~ sample_data(platyH)$Type)

# Hantu Goniopora tissue/skeleton DNA, p=0.04
rich = estimate_richness(goni, measures=c("Shannon", "Simpson", "Chao1"))
t.test(rich$Shannon~ sample_data(goni)$Type)

#-------------------------------
#  RNA Tissue vs Skeleton
#------------------------------

pocil <- subset_samples(ps.all.LR, Species2=="K_Pocil" & NucleicType == "RNA")
platyK <- subset_samples(ps.all.LR, Species2=="K_Platy"& NucleicType == "RNA")
platyH <- subset_samples(ps.all.LR, Species2=="H_Platy"& NucleicType == "RNA")
goni <- subset_samples(ps.all.LR, Species2=="H_Goni"& NucleicType == "RNA")

# Kusu Pocillopora tissue/skeleton RNA, p=0.8
rich = estimate_richness(pocil, measures=c("Shannon", "Simpson"))
t.test(rich$Shannon~ sample_data(pocil)$Type)

# Kusu Platygyra tissue/skeleton RNA, p=0.01
rich = estimate_richness(platyK, measures=c("Shannon", "Simpson"))
t.test(rich$Shannon~ sample_data(platyK)$Type)

# Hantu Platygyra tissue/skeleton RNA, p=0.35
rich = estimate_richness(platyH, measures=c("Shannon", "Simpson"))
t.test(rich$Shannon~ sample_data(platyH)$Type)

# Hantu Goniopora tissue/skeleton RNA, p=0.5
rich = estimate_richness(goni, measures=c("Shannon", "Simpson"))
t.test(rich$Shannon~ sample_data(goni)$Type)



#-------------------------------
#  PERMANOVAS 
#-------------------------------


##---------------------------------------------------------------------
# Tissue vs Skeleton PERMANOVA
##---------------------------------------------------------------------


pocil.all <- subset_samples(ps.all.norm, Species2=="K_Pocil")
platyK.all <- subset_samples(ps.all.norm, Species2=="K_Platy")
platyH.all <- subset_samples(ps.all.norm, Species2=="H_Platy")
goni.all <- subset_samples(ps.all.norm, Species2=="H_Goni")


#---------------------------------
# Pocillopora
#---------------------------------
ps.stats <- pocil.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$SpeciesType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
# Test whether the skeleton and tissue differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:
adonis(BC.dist ~ sample_data(ps.stats)$Type,permutations = 1000)



#---------------------------------
# Platygyra Kusu
#---------------------------------
ps.stats <- platyK.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$SpeciesType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
# Test whether the skeleton and tissue differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:
adonis(BC.dist ~ sample_data(ps.stats)$Type,permutations = 1000)


#---------------------------------
# Platygyra Hantu
#---------------------------------
ps.stats <- platyH.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$SpeciesType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
# Test whether the skeleton and tissue differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:
adonis(BC.dist ~ sample_data(ps.stats)$Type,permutations = 1000)


#---------------------------------
# Goniopora
#---------------------------------
ps.stats <- goni.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix
J.dist = phyloseq::distance(ps.stats, method="jaccard", weighted=F)  #jaccard distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$SpeciesType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
# Test whether the skeleton and tissue differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:
adonis(BC.dist ~ sample_data(ps.stats)$Type,permutations = 1000)


##---------------------------------------------------------------------
# DNA vs RNA PERMANOVA
##---------------------------------------------------------------------


pocil.all <- subset_samples(ps.all.norm, Species2=="K_Pocil")
platyK.all <- subset_samples(ps.all.norm, Species2=="K_Platy")
platyH.all <- subset_samples(ps.all.norm, Species2=="H_Platy")
goni.all <- subset_samples(ps.all.norm, Species2=="H_Goni")

#---------------------------------
# Pocillopora
#---------------------------------
ps.stats <- pocil.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$NucleicType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
adonis(BC.dist ~ sample_data(ps.stats)$NucleicType,permutations = 1000)

#---------------------------------
# Platygyra Kusu
#---------------------------------
ps.stats <- platyK.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$NucleicType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
adonis(BC.dist ~ sample_data(ps.stats)$NucleicType,permutations = 1000)


#---------------------------------
# Platygyra Hantu
#---------------------------------
ps.stats <- platyH.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$NucleicType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
adonis(BC.dist ~ sample_data(ps.stats)$NucleicType,permutations = 1000)


#---------------------------------
# Goniopora
#---------------------------------
ps.stats <- goni.all  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$NucleicType)
permutest(disp.age, pairwise=TRUE, permutations=1000)
adonis(BC.dist ~ sample_data(ps.stats)$NucleicType,permutations = 1000)



##---------------------------------------------------------------------
##----------MAIN MANUSCRIPT FIGURES-------------------------------
##---------------------------------------------------------------------


##---------------------------------------------------------------------
# Figure 3a. Principal Coordinates Analysis
##---------------------------------------------------------------------

logt  = transform_sample_counts(ps.coral.all.norm, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
pcoa.plot1 <-plot_ordination(logt, out.pcoa.logt, type = "samples", shape = "Type2", axes = 1:2, color = "SpeciesType2") + labs(col = "SpeciesType2")+  stat_ellipse(type = "t") +
  scale_color_manual(values=c("#762a83", "#9970ab","#dfc27d", "#bf812d","#92c5de", "#4393c3", "#b8e186", "#7fbc41")) + geom_point(size = 4) + 
  theme_bw() + ggtitle("PCOA coral DNA log transformed by SpeciesType")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10))


# ggsave("all_PCOA_all_ellipse.pdf", plot = pcoa.plot1,  path = "/path",
#        width = 8,
#        height = 5)



##---------------------------------------------------------------------
# Figure 3b. RNA DNA Unifrac Principal Coordinates Analysis
##---------------------------------------------------------------------

pocil.all <- subset_samples(ps.all.norm, Species2=="K_Pocil")
platyK.all <- subset_samples(ps.all.norm, Species2=="K_Platy")
platyH.all <- subset_samples(ps.all.norm, Species2=="H_Platy")
goni.all <- subset_samples(ps.all.norm, Species2=="H_Goni")


logt  = transform_sample_counts(pocil.all, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, "PCoA", "unifrac", weighted=TRUE)
pocil.pcoa1 <- plot_ordination(logt, out.pcoa.logt, color="Type", shape="SpeciesType") +
  stat_ellipse(type = "t")  + scale_color_manual(values=c("#7fbc41","#b8e186","#7fbc41","#b8e186")) + geom_point(size = 4) + scale_shape_manual(values = c(16,1,17,2)) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10))+
  stat_ellipse(type = "t") + ggtitle("Kusu Pocillopora, weighted")


logt  = transform_sample_counts(goni.all, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, "PCoA", "unifrac", weighted=TRUE)
goni.pcoa1 <-plot_ordination(logt, out.pcoa.logt, color="Type", shape="SpeciesType") + stat_ellipse(type = "t")  + scale_color_manual(values=c("#9970ab","#762a83","#9970ab","#762a83")) + geom_point(size = 4) + scale_shape_manual(values = c(16,1,17,2)) +  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10))+  stat_ellipse(type = "t") +
  ggtitle("Hantu Goniopora, weighted")


logt  = transform_sample_counts(platyK.all, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, "PCoA", "unifrac", weighted=TRUE)
platyK.pcoa1 <-plot_ordination(logt, out.pcoa.logt, color="Type", shape="SpeciesType") + stat_ellipse(type = "t")  + scale_color_manual(values=c("#4393c3","#92c5de","#4393c3","#92c5de")) + geom_point(size = 4) +theme_bw() + scale_shape_manual(values = c(16,1,17,2))   + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10))+  stat_ellipse(type = "t") +
  ggtitle("Kusu Platygyra, weighted")


logt  = transform_sample_counts(platyH.all, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, "PCoA", "unifrac", weighted=TRUE)
platyH.pcoa1 <-plot_ordination(logt, out.pcoa.logt, color="Type", shape="SpeciesType") + stat_ellipse(type = "t")  + scale_color_manual(values=c("#bf812d","#dfc27d","#bf812d","#dfc27d")) + geom_point(size = 4) + theme_bw()+ scale_shape_manual(values = c(16,1,17,2)) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10))+  stat_ellipse(type = "t") +
  ggtitle("Hantu Platygyra, weighted")

pcoa.weighted  <- grid.arrange(platyH.pcoa1, goni.pcoa1, pocil.pcoa1, platyK.pcoa1,   nrow = 2)

# 
# ggsave("species_PCOA_unfrac_weighted_updated.pdf", plot = pcoa.weighted,  path = "/path",
#        width = 8,
#        height = 5)


##---------------------------------------------------------------------
# Figure 7. RNA:DNA ratio plot
##---------------------------------------------------------------------


#Make subsets for each species (includind RNA and DNA)
pocil <- subset_samples(ps.coral.all.norm, Species2=="K_Pocil")
platyK <- subset_samples(ps.coral.all.norm, Species2=="K_Platy")
platyH <- subset_samples(ps.coral.all.norm, Species2=="H_Platy")
goni <- subset_samples(ps.coral.all.norm, Species2=="H_Goni")

#Agglomerate to the order level
pocil.glom <- tax_glom(pocil, taxrank="Order", NArm=TRUE)
platyK.glom <- tax_glom(platyK, taxrank="Order", NArm=TRUE)
platyH.glom <- tax_glom(platyH, taxrank="Order", NArm=TRUE)
goni.glom <- tax_glom(goni, taxrank="Order", NArm=TRUE)

#Normalize
pocil.glom.norm <- ps_normalize_median(pocil.glom, "merged_pocil_norm")
platyK.glom.norm <- ps_normalize_median(platyK.glom, "merged_platyK_norm")
platyH.glom.norm <- ps_normalize_median(platyH.glom, "merged_platyH_norm")
goni.glom.norm <- ps_normalize_median(goni.glom, "merged_goniH_norm")

#Merge DNA replicates, merge RNA replicates
pocil.merge <- merge_samples(pocil.glom.norm, "Type")  
platyK.merge <- merge_samples(platyK.glom.norm, "Type")  
platyH.merge <- merge_samples(platyH.glom.norm, "Type")  
goni.merge <- merge_samples(goni.glom.norm, "Type")  

#Normalize again
pocil.merge.norm <- ps_normalize_median(pocil.merge, "merged_pocil_norm")  
platyK.merge.norm <- ps_normalize_median(platyK.merge, "merged_platyK_norm")  
platyH.merge.norm <- ps_normalize_median(platyH.merge, "merged_platyH_norm")  
goni.merge.norm <- ps_normalize_median(goni.merge, "merged_goniH_norm")  


#Convert otu table to dataframe
pocil.df <-phyloseq_to_df(pocil.merge.norm, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance")
platyK.df <-phyloseq_to_df(platyK.merge.norm, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance")
platyH.df <-phyloseq_to_df(platyH.merge.norm, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance")
goni.df <-phyloseq_to_df(goni.merge.norm, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance")


#Calculate RNA/DNA ratios
pocil.df <- transform(pocil.df, T_ratio = Tissue_R /Tissue, S_ratio=Skeleton_R/Skeleton)
platyK.df <- transform(platyK.df, T_ratio = Tissue_R /Tissue, S_ratio=Skeleton_R/Skeleton)
platyH.df <- transform(platyH.df, T_ratio = Tissue_R /Tissue, S_ratio=Skeleton_R/Skeleton)
goni.df <- transform(goni.df, T_ratio = Tissue_R /Tissue, S_ratio=Skeleton_R/Skeleton)


#Reorganize dataframe to faclitate plot_bar
#----------------------------------------------------------------------------------------------------------------------
#------Pocillopora-------#
pocil.df.S <- as_tibble(pocil.df) %>% select(Phylum, Class, Order, Skeleton, Skeleton_R, S_ratio)%>% transform(., Type="Skeleton", Coral="Kusu Pocillopora")%>% dplyr::rename(DNA=Skeleton, RNA=Skeleton_R, Ratio=S_ratio)
pocil.df.T <- pocil.df %>% select(Phylum, Class, Order,Tissue, Tissue_R, T_ratio)%>% transform(., Type="Tissue", Coral="Kusu Pocillopora")%>% dplyr::rename(DNA=Tissue, RNA=Tissue_R, Ratio=T_ratio)
pocil.df2 <-as_tibble(rbind(pocil.df.T, pocil.df.S))%>%na.omit()
#----------------------------------------------------------------------------------------------------------------------
#------Platygyra Kusu-------#
platyK.df.S <- as_tibble(platyK.df) %>% select(Phylum, Class, Order, Skeleton, Skeleton_R, S_ratio)%>% transform(., Type="Skeleton", Coral="Kusu Platgyra")%>% dplyr::rename(DNA=Skeleton, RNA=Skeleton_R, Ratio=S_ratio)
platyK.df.T <- platyK.df %>% select(Phylum, Class, Order,Tissue, Tissue_R, T_ratio)%>% transform(., Type="Tissue", Coral="Kusu Platgyra")%>% dplyr::rename(DNA=Tissue, RNA=Tissue_R, Ratio=T_ratio)
platyK.df2 <-as_tibble(rbind(platyK.df.T, platyK.df.S))%>%na.omit()
#----------------------------------------------------------------------------------------------------------------------
#------Platygyra Hantu-------#
platyH.df.S <- as_tibble(platyH.df) %>% select(Phylum, Class, Order, Skeleton, Skeleton_R, S_ratio)%>% transform(., Type="Skeleton", Coral="Hantu Platgyra")%>% dplyr::rename(DNA=Skeleton, RNA=Skeleton_R, Ratio=S_ratio)
platyH.df.T <- platyH.df %>% select(Phylum, Class, Order,Tissue, Tissue_R, T_ratio)%>% transform(., Type="Tissue", Coral="Hantu Platgyra")%>% dplyr::rename(DNA=Tissue, RNA=Tissue_R, Ratio=T_ratio)
platyH.df2 <-as_tibble(rbind(platyH.df.T, platyH.df.S))%>%na.omit()
#----------------------------------------------------------------------------------------------------------------------
#------Goniopora-------#
goni.df.S <- as_tibble(goni.df) %>% select(Phylum, Class, Order, Skeleton, Skeleton_R, S_ratio)%>% transform(., Type="Skeleton", Coral="Hantu Goniopora")%>% dplyr::rename(DNA=Skeleton, RNA=Skeleton_R, Ratio=S_ratio)
goni.df.T <- goni.df %>% select(Phylum, Class, Order,Tissue, Tissue_R, T_ratio)%>% transform(., Type="Tissue", Coral="Hantu Goniopora")%>% dplyr::rename(DNA=Tissue, RNA=Tissue_R, Ratio=T_ratio)
goni.df2 <-as_tibble(rbind(goni.df.T, goni.df.S))%>%na.omit()
#----------------------------------------------------------------------------------------------------------------------

#Combine dataframes 
all.coral.ratios <- as_tibble(rbind(goni.df2, platyH.df2,platyK.df2,pocil.df2))

#----------------------------------------------------------------------------------------------------------------------------

#Select orders that were identified in nifH analysis
select.ratios <- filter(all.coral.ratios, Order == "Synechococcales" | Order == "Oscillatoriales" | Order == "Chroococcales" | Order =="Pleurocapsales" | Order =="Chrococcales" | Order =="Nostocales" |
                          Order =="Caldilineales" | Order =="Ardenticatenales" | Order == "Anaerolineales"| Order =="Chloroflexales" | Order =="SAR202 clade" |  Order =="SBR1031" |
                          Order =="Oceanospirillales" | Order =="Alteromonadales" | Order =="Thiotrichales" | Order == "Vibrionales" | Order =="Nitrosococcales"| Order =="Chromatiales" | 
                          Order =="Kiloniellales" | Order =="Rhodospirillales" | Order == "Rhizobiales" | Order == "Rhodobacterales" | Order== "Desulfobacterales" | Order == "Desulfobulbales" | 
                          Order =="Desulfovibrionales" | Order == "Syntrophobacterales" |  Class == "Nitrospiria" )

# Save as CSV files 
# write.csv(select.ratios, file = "/Users/molly/Dropbox/NitrogenFixation_Singapore/Manuscripts/August_experiment_manuscript/overleaf_supplemental/select_RNA_DNA_ratios_byspecies.csv")
# select.ratios.test <- filter(all.coral.ratios, Order == "Synechococcales" | Order == "Oscillatoriales" | Order == "Chroococcales" | Order =="Pleurocapsales" )

#Make new column with taxonomic ranks combined, so that samples will be sorted in a logical order on the y-axis
select.ratios2 <- unite(select.ratios, "Phylum_Class_Order", Phylum:Class:Order,sep = "_",  remove = FALSE )

finalplot_log <- ggplot(select.ratios2, aes(x= Ratio, y= factor(Phylum_Class_Order, levels=rev(levels(factor(Phylum_Class_Order)))), shape=Coral, fill = factor(Type)))+
  geom_point(size = 3) +  scale_shape_manual(values=c(22, 24, 25,21))+
  scale_fill_manual(values=c("#fee08b", "#8073ac"))+
  theme_bw() +  ggtitle("RNA:DNA")+ xlab("log RNA:DNA ratio") + ylab("Order") +
  scale_x_continuous(trans = 'log10') +
  annotation_logticks(sides="b") 

# ggsave("select_RNA_DNA_ratios_byspecies_log.pdf", plot = finalplot_log, path = "/path",
#        width = 10,
#        height = 9)



##---------------------------------------------------------------------
##----------SUPPLEMENTAL FIGURES-------------------------------
##---------------------------------------------------------------------




##---------------------------------------------------------------------
# Figre S6. Treemaps 
##---------------------------------------------------------------------


# Deine function to convert phyloseq for treemaps
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




# Function fortTreemap with custom color palette 
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


# Make custom color palette for the order level 
getPalette = colorRampPalette(c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))) 
orderList = unique(tax_table(ps.coral.all.abund)[,"Order"])
orderPalette = getPalette(length(orderList))
names(orderPalette) = orderList

#Make subset of data for different teemaps
tree.dna <- subset_samples(ps.coral.all.abund, NucleicType=="DNA")
tree.rna <- subset_samples(ps.coral.all.abund, NucleicType=="RNA")
tissue.dna <- subset_samples(tree.dna, Type2=="Tissue")
skeleton.dna <- subset_samples(tree.dna, Type2=="Skeleton")
tissue.rna <- subset_samples(tree.rna, Type2=="Tissue")
skeleton.rna <- subset_samples(tree.rna, Type2=="Skeleton")


#convert to ps_to_long 
long_tissue_rna <- ps_to_long(tissue.rna)
long_skeleton_rna <- ps_to_long(skeleton.rna)
long_tissue_dna <- ps_to_long(tissue.dna)
long_skeleton_dna <- ps_to_long(skeleton.dna)

#DNA tissue treemap
dna.tissue.tree.order <- treemap_gg_dv2(long_tissue_dna, Class, Order,"Coral 16S rRNA community - Tissue DNA samples", orderPalette)
#DNA skeleton treemap
dna.skeleton.tree.order <- treemap_gg_dv2(long_skeleton_dna, Class, Order,"Coral 16S rRNA community - Skeleton DNA samples",orderPalette)
#RNA tissue treemap
rna.tissue.tree.order <- treemap_gg_dv2(long_tissue_rna, Class, Order,"Coral 16S rRNA community - Tissue RNA samples", orderPalette)
#RNA skeleton treemap
rna.skeleton.tree.order <- treemap_gg_dv2(long_skeleton_rna, Class, Order,"Coral 16S rRNA community - Skeleton RNA samples",orderPalette)

x.order <-grid.arrange(dna.tissue.tree.order, dna.skeleton.tree.order, rna.tissue.tree.order, rna.skeleton.tree.order, nrow = 2)

# ggsave("ST_DNA_RNA_treeplots_order_updated.pdf", plot = x.order,
#        path = "/path",width = 10, height = 8)



##---------------------------------------------------------------------
# Figure S7. Bar plots 
##---------------------------------------------------------------------

#Fig S7a
top <- names(sort(taxa_sums(ps.coral.dna.abund), decreasing=TRUE))[1:300]
ps.top <- transform_sample_counts(ps.coral.dna.abund, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(top, ps.top)
dna.bar <- plot_bar(ps.top, fill="Order", x="Sample") + facet_wrap(~SpeciesType, 2, scales="free_x")+
  theme(legend.key.size = unit(0.9,"line"))  + theme(legend.position="bottom")

#Fig S7b
top <- names(sort(taxa_sums(ps.coral.rna.abund), decreasing=TRUE))[1:300]
ps.top <- transform_sample_counts(ps.coral.rna.abund, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(top, ps.top)
rna.bar <- plot_bar(ps.top, fill="Order", x="Sample") + facet_wrap(~SpeciesType, 2, scales="free_x")+
  theme(legend.key.size = unit(0.9,"line"))  + theme(legend.position="bottom")


##---------------------------------------------------------------------
# Figure S12. Seawater community 16S tree map
##---------------------------------------------------------------------

# Define color palettes for Seawater
getPalette = colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))) 
orderList = unique(tax_table(ps.SW.norm)[,"Order"])
orderPalette = getPalette(length(orderList))
names(orderPalette) = orderList


ps_to_long_SW <- function(ps) {
  otu_df <- data.frame(otu_table(ps)) %>% 
    rownames_to_column(var = "otu_id")
  taxo_df <- data.frame(tax_table(ps))%>% 
    rownames_to_column(var = "otu_id")
  otu_df <- left_join(taxo_df, otu_df)
  otu_df <- gather(otu_df, "sample","n_seq" ,contains("SW")) # All samples contain X
  metadata_df <- data.frame(sample_data(ps)) %>% 
    rownames_to_column(var = "sample")
  otu_df <- left_join(otu_df, metadata_df)
}

# Seawater tree maps
long_SW <- ps_to_long_SW(ps.SW.norm)
dna.tree.order <- treemap_gg_dv2(long_SW, Class, Order,"Seawater 16S DNA-based community", orderPalette)

# ggsave("SW_16S_community.pdf", plot = dna.tree.order,
#        path = "/path",
#        width = 5,
#        height = 5)





