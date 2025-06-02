#
rm(list = ls())
gc()
suppressMessages({
  library("tidyverse")
  library("phyloseq")
  library("rstatix")
  library("vegan")
  library("picante")
  library("ggpubr")
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(poppr)  # For AMOVA if genetic analysis is needed
  library(microbiomeutilities)
  library(viridis)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(ggtext)
})
no_of_cores = 16
setwd("/work/benson/bpeng4/Kemin/Kemin_All")
load("intermediate/Kemin_ps.rda")
ls()
ps

###Pre-arrangement on  the phyloseq file
################
ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))
ps.clean
#Filter Out Taxa Exist in At Least 25% samples
ps.clean.p0 <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 282}, prune=TRUE)
ps.clean.p0

#Complete Metadata
sample_meta = sample_data(ps.clean.p0)
sample_meta$Rep <- as.character(sample_meta$Rep)
sample_meta$MicroSam<- paste(sample_meta$Microbiome,sample_meta$Sample, sep = "_")
sample_meta$MicroCatg<- paste(sample_meta$Microbiome,sample_meta$Category., sep = "_")
sam_data(ps.clean.p0)<-sample_meta
sam_data(ps)<-sample_meta

#Check Total Reads to find the Thresholds for Rarafication
OTU<-otu_table(ps.clean.p0) |>as.data.frame()
OTU$Total<- rowSums(OTU)

#Rarefication
ps.rare<-rarefy_even_depth(ps.clean.p0, sample.size = 12487,
                           rngseed = 111, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
#Change to relative abundance
ps.clean.re <- transform_sample_counts(ps.rare, function(x) x / sum(x))
ps.sample.re<-subset_samples(ps.clean.re, !ps.clean.re@sam_data$Sample %in% c("FBB0","FBB16"))
######################

###Plot Stack Bar chart for Baseline relative abundance
####################
ps.re.FBB<-subset_samples(ps.clean.re, ps.clean.re@sam_data$Sample %in% c("FBB0","FBB16"))

#Set up mycolor
mycolor<-c(
  '#F781BF', '#CCEBC5', '#FFFF99', '#FF7F00', '#8DD3C7', '#7FC97F', '#999999',  
  '#CAB2D6', '#FDBF6F', '#D95F02', '#666666', '#A6761D', '#B2DF8A',  
  '#6A3D9A', '#FB9A99', '#E41A1C', '#FFFFB3', '#FDB462', '#F0027F', '#D9D9D9', '#4DAF4A', 
  '#FFFF33', '#B15928', '#FDC086', '#80B1D3', '#A6CEE3', '#BC80BD', '#1B9E77', '#E31A1C',  
  '#FCCDE5', '#386CB0', '#377EB8', '#984EA3', '#FFED6F', '#66A61E', '#E6AB02',  
  '#BF5B17', '#1F78B4', '#A65628', '#B3DE69', '#7570B3', '#E7298A', '#33A02C','#FB8072',
  '#f7754f', '#ce9032', '#97a431', '#32b166', '#35ad9c', '#38a9c5', '#3ca2f4','#a48cf4', 
  '#f45cf2', '#f66bad')

#Plot for Microbiome Baseline Profile
phyloseq::plot_bar(ps.re.FBB, fill = "Family") + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~MicroSam, scales = "free_x",ncol = 4) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, size = 8),
        legend.key.size = unit(17, "pt"))+
  scale_fill_manual(values=mycolor)+
  guides(fill = guide_legend(ncol = 1))

########################

####Beta Diversity Analysis and Plotting
########################
ps.beta= ps.clean.re

###Principle Component Analysis
ordBC <- ordinate(ps.beta, "PCoA", "bray")
ordJC <- ordinate(ps.beta, "PCoA", "jaccard")
ordUF <- ordinate(ps.beta, "PCoA", "unifrac")
ordwUF <- ordinate(ps.beta, "PCoA", "wunifrac")
smpID <- sample_data(ps.beta)$MicroSam

# add sample_data info
BC <- merge(data.frame(ordBC$vectors[,1:2], sample = smpID, method = 'BC'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(BC) <- BC$Row.names
BC$Row.names <- NULL  # Remove the extra column created by merge

Jaccard <- merge(data.frame(ordJC$vectors[,1:2], sample = smpID,method = 'Jaccard'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(Jaccard) <- Jaccard$Row.names
Jaccard$Row.names <- NULL  # Remove the extra column created by merge

unifrac <- merge(data.frame(ordUF$vectors[,1:2], sample = smpID,method = 'unifrac'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(unifrac) <- unifrac$Row.names
unifrac$Row.names <- NULL  # Remove the extra column created by merge

wunifrac <- merge(data.frame(ordwUF$vectors[,1:2], sample = smpID,method = 'wunifrac'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(wunifrac) <- wunifrac$Row.names
wunifrac$Row.names <- NULL  # Remove the extra column created by merge

# Keep first 2 vectors (latent variables, PCs) of each distance matrix
df <- rbind(BC,Jaccard,unifrac,wunifrac)

# Calculate the mean for each level of MicroSample
summary_samdf <- df %>%
  group_by(method,Microbiome,Sample, MicroSam) %>%
  summarise(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE)
  )
#Order Sample levels
summary_samdf$Sample <- factor(summary_samdf$Sample, 
                               levels = c("BetaVia Pure WD", "Oat Bagasse hydrolysate-(B)",         
                                          "Oat Bagasse-Control", "Oat Bagasse hydrolysate-(PBV)",       
                                          "beetroot powder", "carrot powder",                       
                                          "dried spinich powder", "ginger powder",                       
                                          "lemon powder", "Olive Flour",                         
                                          "pomegranate powder", "Cellulose",                           
                                          "Ground Beet Pulp", "HiSmooth flax seed fiber",            
                                          "Micronized cassava fiber", "Psyllium Husk",                       
                                          "guar gum", "kappa carageenan ",                   
                                          "locust bean gum", "xanthan gum",                         
                                          "Jerusalum articohoke powder (inulin)", "Mushroom chitosan",                   
                                          "dried apple pomace","FBB16","FBB0"))
#Set up FBB0 as the starting points
summary_samdf0<-summary_samdf
summary_samdf0<-summary_samdf0[summary_samdf0$Sample!="FBB16",]
fbb0_samdf <- summary_samdf0 %>%
  filter(Sample == "FBB0") %>%
  select(method, Microbiome, Axis.1, Axis.2) %>%
  mutate(Axis.1_start = Axis.1, Axis.2_start = Axis.2) %>%
  select(-Axis.1, -Axis.2)  # Remove old columns to avoid duplicates
summary_samdf0 <- summary_samdf0 %>%
  left_join(fbb0_samdf, by = c("Microbiome", "method"))  # Match on method & Microbiome

#Plot with arrows
ggplot(summary_samdf0, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal()


#Beta Diversity Plot with arrows for each microbiome
####################
summary_df0_K3<-summary_samdf0[summary_samdf0$Microbiome=="K3",]
summary_df0_K5<-summary_samdf0[summary_samdf0$Microbiome=="K5",]
summary_df0_K7<-summary_samdf0[summary_samdf0$Microbiome=="K7",]
summary_df0_K10<-summary_samdf0[summary_samdf0$Microbiome=="K10",]
#Get the first four colors
first_four <- brewer.pal(n = 8, name = "Set2")[1:4]
ggplot(summary_df0_K10, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[1]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K10", "Beta Diversity Shift from the Baseline"))
ggplot(summary_df0_K3, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[2]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K3", "Beta Diversity Shift from the Baseline"))
ggplot(summary_df0_K5, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[3]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K5", "Beta Diversity Shift from the Baseline"))
ggplot(summary_df0_K7, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[4]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K7", "Beta Diversity Shift from the Baseline"))
#############

#Beta Diversity Plot with arrows for each group of health and disease microbiomes
####################
summary_df0<-summary_samdf0
summary_df0_K3K5<-summary_df0[summary_df0$Microbiome%in%c("K3","K5"),]
summary_df0_K3K7<-summary_df0[summary_df0$Microbiome%in%c("K3","K7"),]
summary_df0_K10K3<-summary_df0[summary_df0$Microbiome%in%c("K10","K3"),]
summary_df0_K10K5<-summary_df0[summary_df0$Microbiome%in%c("K10","K5"),]
summary_df0_K10K7<-summary_df0[summary_df0$Microbiome%in%c("K10","K7"),]
#Get the first four colors
first_four <- brewer.pal(n = 8, name = "Set2")[1:4]
ggplot(summary_df0_K3K5, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(2,3)]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K3&K5", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K3K7, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(2,4)]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K3&K7", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K10K3, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(1,2)]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K10&K3", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K10K5, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(1,3)]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K10&K5", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K10K7, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(1,4)]) +
  scale_shape_manual(values = c(0:18,20:23,19)) +  # Adjust shapes as needed
  labs(shape = "Sample", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K10&K7", "Beta Diversity Shift from the Baseline"))

