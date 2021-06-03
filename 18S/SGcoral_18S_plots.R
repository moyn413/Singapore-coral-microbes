#----------------------------------------------------------------
# 18S phyloseq processing file and script to make plots with 18S data. Starting file is a phyloseq file made using the dada2 pipeline.
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


ps <- readRDS("/GitHub/Singapore-coral-microbes/18S/SGcoral_18S_phyloseq.rds")   #read RDS phyloseq file

#--------------------------------
# Remove contaminants    
#--------------------------------

#Using decontam function, identify contaminants using the prevalence method with a threshold of 0.1
sample_data(ps)$is.neg <- sample_data(ps)$Nucleic == "detcontam"
contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant) #no contaminants identified
# remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev01$contaminant, ps) 

#remove samples used for decontam
ps.noncontam_noC <-  ps.noncontam %>% subset_samples(Nucleic != "decontam") 

#remove metazoans
ps.noMet <- subset_taxa(ps.noncontam_noC, Division != "Metazoa"| is.na(Division)) 

# Make separate phyloseqs for different sample types
ps.all <- ps.noMet
ps.coral.all <-  ps.noMet %>% subset_samples(Type != "Seawater")
ps.coral.dna <-  ps.noMet %>% subset_samples(Type != "Seawater") %>% subset_samples(Nucleic == "DNA")
ps.coral.rna <-  ps.noMet %>% subset_samples(Type != "Seawater") %>% subset_samples(Nucleic == "RNA")
ps.SW <-  ps.noMet %>% subset_samples(Type == "Seawater")


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
ps.SW.LO <- prune_taxa(taxa_sums(ps.SW.LR)>=10, ps.SW.LR)  


#----------------------------------------------------------------------------------------------
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
#----------------------------------------------------------------------------------------------

#--------------------------------
# Normalize phyloseqs
#--------------------------------

# Normalize phyloseq files with low abundant OTUs (for diversity plots)
ps.all.LR.norm = ps_normalize_median(ps.all.LR, "18S_all_LR_norm")
ps.coral.all.LR.norm = ps_normalize_median(ps.coral.all.LR, "Coral_18S_all_LR_norm")
ps.coral.dna.LR.norm = ps_normalize_median(ps.coral.dna.LR, "Coral_18S_dna_LR_norm")
ps.coral.rna.LR.norm = ps_normalize_median(ps.coral.rna.LR, "Coral_18S_rna_LR_norm")
ps.SW.LR.norm = ps_normalize_median(ps.SW.LR, "SW_18S_norm")

# Normalize phyloseq file without low abundant OTUs
ps.all.norm = ps_normalize_median(ps.all.LO, "18S_all_norm")
ps.coral.all.norm = ps_normalize_median(ps.coral.all.LO, "Coral_18S_all_norm")
ps.coral.dna.norm = ps_normalize_median(ps.coral.dna.LO, "Coral_18S_dna_norm")
ps.coral.rna.norm = ps_normalize_median(ps.coral.rna.LO, "Coral_18S_rna_norm")
ps.SW.norm = ps_normalize_median(ps.SW.LO, "SW_18S_norm")

# Remove non-abundant taxa and normalize again
ps.all.abund = ps_abundant(ps.all.norm,contrib_min=0.01, "18S_all_abund")  
ps.coral.all.abund = ps_abundant(ps.coral.all.norm,contrib_min=0.01, "Coral_18S_all_abund")  
ps.coral.dna.abund = ps_abundant(ps.coral.dna.norm,contrib_min=0.01, "Coral_18S_dna_abund")  
ps.coral.rna.abund = ps_abundant(ps.coral.rna.norm,contrib_min=0.01, "Coral_18S_rna_abund")  
ps.SW.abund = ps_abundant(ps.SW.norm,contrib_min=0.003, "SW_18S_abund") 


# Remove Symbiodinium and noramlize again 
ps.all.nosym <- subset_taxa(ps.all.LO, Genus != "Symbiodinium"| is.na(Genus)) %>% ps_normalize_median(.,"18S_all_norm" )
ps.coral.nosym <- subset_taxa(ps.coral.all.LO, Genus != "Symbiodinium"| is.na(Genus)) %>%ps_normalize_median(.,"Coral_18S_all_norm" )
ps.coral.dna.nosym <- subset_taxa(ps.coral.dna.LO, Genus != "Symbiodinium"| is.na(Genus)) %>% ps_normalize_median(.,"Coral_18S_dna_norm" )
ps.coral.rna.nosym <- subset_taxa(ps.coral.rna.LO, Genus != "Symbiodinium"| is.na(Genus)) %>% ps_normalize_median(.,"Coral_18S_rna_norm" )

ps.coral.nosym.LO <- subset_taxa(ps.coral.all.LO, Genus != "Symbiodinium"| is.na(Genus))



#--------------------------------
# Alpha diversity stats
#--------------------------------

# Perform stats on all samples, not normalized, not trimmed. Use this for stats table
physeq2 = transform_sample_counts(ps.all.LR, ceiling)
rich = estimate_richness(physeq2)
plot_richness(ps.all, measures=c("Shannon", "Simpson"), x="SpeciesType", color="Type") + geom_boxplot()

# DNA tissue vs skeleton 
plot_richness(ps.coral.dna.LR.norm, measures=c("Shannon", "Simpson"), x="SpeciesType", color="Type") + geom_boxplot()
rich = estimate_richness(ps.coral.dna.LR.norm, measures=c("Shannon", "Simpson"))
# Shannon stats, Wilcoxon rank sum test and t.test
t.test(rich$Shannon~ sample_data(ps.coral.dna.LR.norm)$Type)


# RNA tissue vs skeleton 
plot_richness(ps.coral.rna.LR.norm, measures=c("Shannon", "Simpson"), x="SpeciesType2", color="Type2") + geom_boxplot()
rich = estimate_richness(ps.coral.rna.LR.norm, measures=c("Shannon", "Simpson"))
# Shannon stats, Wilcoxon rank sum test and t.test
t.test(rich$Shannon~ sample_data(ps.coral.rna.LR.norm)$Type)


#--------------------------------
# PERMANOVA
#--------------------------------

ps.stats <- ps.coral.all.LO  # perform stats on normalized file, with OTUs < 10x removed
BC.dist = phyloseq::distance(ps.stats, method="bray", weighted=F)  #bray distance matrix

# Beta dispersion- first test that dispersions are the same. If all results are non-significant, proceed with PERMANOVA 
disp.age = betadisper(BC.dist, sample_data(ps.stats)$Species2)
permutest(disp.age, pairwise=TRUE, permutations=1000)
adonis(BC.dist ~ sample_data(ps.stats)$Species,permutations = 1000)



##---------------------------------------------------------------------
##----------Supplemental Figures-------------------------------
##---------------------------------------------------------------------


##---------------------------------------------------------------------
# Figure S8. Symbiodinium relative abundance
##---------------------------------------------------------------------

getPalette = colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(8, "Set1"))) 
genusList = unique(tax_table(ps.coral.all.abund)[,"Species"])
genusPalette = getPalette(length(genusList))
names(genusPalette) = genusList


# Symbiodinium relative abundance bar plot by genus
coral.tissue <- subset_samples(ps.coral.all.abund, Type2 == "Tissue")
top_c <- names(sort(taxa_sums(coral.tissue), decreasing=TRUE))[1:200]
ps.top_c <- transform_sample_counts(coral.tissue, function(OTU) OTU/sum(OTU))
ps.top_c <- prune_taxa(top_c , ps.top_c)
genus <- plot_bar(ps.top_c, fill="Species", x="Name") + facet_wrap(~Type, 2, scales="free_x")+
  theme(legend.key.size = unit(0.8,"line"))  + theme(legend.position="bottom") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + scale_color_manual(values=genusPalette) + scale_fill_manual(values= genusPalette)


# ggsave("18S_coral.pdf", plot = genus, path = "/path",
#        width = 10,
#        height = 7)




##---------------------------------------------------------------------
# Figure S9. 18S Principal coordinates analysis
##---------------------------------------------------------------------

logt  = transform_sample_counts(ps.coral.all.norm, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
pcoa.plot1 <-plot_ordination(logt, out.pcoa.logt, type = "samples", shape = "Type2", axes = 1:2, color = "SpeciesType2") + labs(col = "SpeciesType2") +
  #  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  scale_color_manual(values=c("#762a83", "#9970ab","#dfc27d", "#bf812d","#92c5de", "#4393c3", "#b8e186", "#7fbc41")) + geom_point(size = 4) + 
  theme_bw() + ggtitle("PCOA coral DNA & RNA log transformed by SpeciesType")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10))
pcoa.plot1

# ggsave("18S_pcoa.pdf", plot = pcoa.plot1,  path = "/path",
#        width = 8,
#        height = 5)


##---------------------------------------------------------------------
# Figure S10. 18S tree maps with no Symbiodinium
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
    theme(legend.position="bottom", plot.title = element_text(size = 16, face = "bold")) #Change legend position to "bottom" for legend
  print(g_treemap)
  return(g_treemap)
}


# Make custom color palette for the order level, * + 12 colors repeat
getPalette = colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))) 
orderList = unique(tax_table(ps.coral.nosym)[,"Order"])
orderPalette = getPalette(length(orderList))
names(orderPalette) = orderList


#Convert phyloseq to long format
tissue.rna <- subset_samples(ps.coral.rna.nosym, Type2=="Tissue")
skeleton.rna <- subset_samples(ps.coral.rna.nosym, Type2=="Skeleton")

tissue.dna <- subset_samples(ps.coral.dna.nosym, Type2=="Tissue")
skeleton.dna <- subset_samples(ps.coral.dna.nosym, Type2=="Skeleton")


#tissue/skeleton
long_tissue_rna <- ps_to_long(tissue.rna)
long_skeleton_rna <- ps_to_long(skeleton.rna)

long_tissue_dna <- ps_to_long(tissue.dna)
long_skeleton_dna <- ps_to_long(skeleton.dna)



# Compare tissue/skeleton DNA
dna.tissue.tree.order <- treemap_gg_dv2(long_tissue_dna, Class, Order ,"Coral 18S rRNA community -  Tissue DNA samples", orderPalette)
dna.skeleton.tree.order <- treemap_gg_dv2(long_skeleton_dna, Class, Order,"Coral 18S rRNA community - Skeleton DNA samples", orderPalette)


# Compare tissue/skeleton RNA
rna.tissue.tree.order <- treemap_gg_dv2(long_tissue_rna, Class, Order ,"Coral 18S rRNA community -  Tissue RNA samples", orderPalette)
rna.skeleton.tree.order <- treemap_gg_dv2(long_skeleton_rna, Class, Order,"Coral 18S rRNA community - Skeleton RNA samples", orderPalette)


# Combine DNA and RNA community plots
x.order <-grid.arrange(dna.tissue.tree.order, dna.skeleton.tree.order, rna.tissue.tree.order, rna.skeleton.tree.order, nrow = 2)


# ggsave("18S_nosymb_treeplots_order.pdf", plot = x.order,
#        path = "/path",width = 10, height = 8)




##---------------------------------------------------------------------
# Figure S12b. Seawater 18S treemaps
##---------------------------------------------------------------------


# Make custom color palette for the family level
getPalette = colorRampPalette(c(brewer.pal(12, "Paired"), brewer.pal(8, "Set1"),brewer.pal(8, "Dark2") )) 
familyList = unique(tax_table(ps.SW.abund)[,"Family"])
familyPalette = getPalette(length(familyList))
names(familyPalette) =familyList


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
#----------------------------#----------------------------
# Seawater tree maps
#----------------------------#----------------------------
long_SW <- ps_to_long_SW(ps.SW.abund)
dna.tree.order <- treemap_gg_dv2(long_SW, Class, Family,"Seawater 18S DNA-based community", familyPalette)
# 
# ggsave("SW_18S_community.pdf", plot = dna.tree.order,
#        path = "/path",
#        width = 5,
#        height = 5)


