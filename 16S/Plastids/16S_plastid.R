library(plyr) 
library(dplyr)
library(DESeq2)
library(ggplot2)
library(vegan)
library(phyloseq)
library(ape)
library(decontam)
library(stringr)
library(tibble)
library(phangorn)
library(seqinr)
library(DECIPHER)
library(dada2)
library(kableExtra)
library(gridExtra)
library(tidyverse)
library(knitr)
library(RColorBrewer)


#--------------------------------
# Load phyloseq file
#--------------------------------


ps <- readRDS("/Users/molly/Dropbox/NitrogenFixation_Singapore/Molecular/SEQUENCING_RESULTS2020/16S/16S_August_coral_SW_plus_wneg_new_silva/dada2/16S_chloroplast_plastid_PR2database/CORAL16S_phyloseq_chloroplast_90.rds")   #read RDS phyloseq file

#--------------------------------
# Remove contaminants    
#--------------------------------

# Identify contaminants using 0.2 threshold and prevalence method
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant)

# remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev01$contaminant, ps) 


# remove control files
ps.noncontam_noeuk_noC <-  ps.noncontam %>% subset_samples(Sample_or_Control != "Control") %>% subset_samples(Sample_or_Control != "noRT") %>% subset_samples(NucleicType != "decontam") 


# Remove OTUs occurring fewer than 10x
ps.all.LO = prune_taxa(taxa_sums(ps.noncontam_noeuk_noC)>=10, ps.noncontam_noeuk_noC)


ps.all <- ps.all.LO
ps.coral.all <-  ps.all.LO %>% subset_samples(Type != "Seawater")
ps.SW.all <-  ps.all.LO %>% subset_samples(Type == "Seawater")



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


#Normalize 
ps.all.norm = ps_normalize_median(ps.all, "16S_all_LR_norm")
ps.coral.all.norm = ps_normalize_median(ps.coral.all, "Coral_16S_all_LR_norm")
ps.SW.norm = ps_normalize_median(ps.SW.all, "Coral_16S_all_LR_norm")









##---------------------------------------------------------------------
# Figure S8. 16S coral plastid sequences
##---------------------------------------------------------------------

getPalette = colorRampPalette(c( brewer.pal(8,"Spectral"), brewer.pal(8,"Dark2"),  brewer.pal(9, "Set1"),brewer.pal(12, "Set3")))
orderList = unique(tax_table(ps.coral.all.norm)[,"Order"])
orderPalette = getPalette(length(orderList))
names(orderPalette) = orderList


# Normalized
top_c <- names(sort(taxa_sums(ps.coral.all.norm), decreasing=TRUE))
ps.top_c <- transform_sample_counts(ps.coral.all.norm, function(OTU) OTU/sum(OTU))
ps.top_c <- prune_taxa(top_c , ps.top_c)
order <-plot_bar(ps.top_c, fill="Order", x="Bars2") + facet_wrap(~Type, 2, scales="free_x")+
  theme(legend.key.size = unit(0.8,"line"))  + theme(legend.position="bottom") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  + scale_color_manual(values=orderPalette) + scale_fill_manual(values= orderPalette)


getPalette = colorRampPalette(c(brewer.pal(9, "Set1"),brewer.pal(12, "Set3"),brewer.pal(8,"Spectral")))
genusList = unique(tax_table(ps.coral.all.norm)[,"Genus"])
genusPalette = getPalette(length(genusList))
names(genusPalette) = genusList


# Abundance, pruned to most abundant taxa
top_c <- names(sort(taxa_sums(ps.coral.all), decreasing=TRUE))[1:50]
ps.top_c <- prune_taxa(top_c , ps.coral.all)
genus <- plot_bar(ps.top_c, fill="Genus", x="Bars2") + facet_wrap(~Type, 2, scales="free_x")+
  theme(legend.key.size = unit(0.8,"line"))  + theme(legend.position="bottom") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + scale_color_manual(values=genusPalette) + scale_fill_manual(values= genusPalette)


x <-grid.arrange(order, genus, nrow = 1)

# ggsave("plastid_barplots.pdf", plot = x, path = "/Users/molly/Dropbox/NitrogenFixation_Singapore/Manuscripts/August_experiment_manuscript/overleaf_supplemental",
#        width = 20,
#        height = 15)



##---------------------------------------------------------------------
# Figure S13a. 16S seawater plastid sequences
##---------------------------------------------------------------------

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

# Define color palettes for Seawater
getPalette = colorRampPalette(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))) 
familyList = unique(tax_table(ps.SW.norm)[,"Family"])
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

# Seawater tree map
long_SW <- ps_to_long_SW(ps.SW.norm)
dna.tree.order <- treemap_gg_dv2(long_SW, Class, Family,"Seawater Plastid 16S DNA-based community", familyPalette)

# ggsave("SW_16S_Plastid_community.pdf", plot = dna.tree.order,
#        path = "/Users/molly/Dropbox/NitrogenFixation_Singapore/Manuscripts/August_experiment_manuscript/overleaf_supplemental",
#        width = 5,
#        height = 5)
