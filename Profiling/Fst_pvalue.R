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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(poppr)  # For AMOVA if genetic analysis is needed
})
no_of_cores = 36
setwd("/work/benson/bpeng4/Kemin/Kemin_All")
load("intermediate/Kemin_ps.rda")
ls()
ps

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

#Plot Stack Bar chart for Baseline relative abundance
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

#####Beta Diversity Analysis and Plotting
########################
ps.beta= ps.clean.re

###Principle Component Analysis
ordBC <- ordinate(ps.beta, "PCoA", "bray")
ordJC <- ordinate(ps.beta, "PCoA", "jaccard")
ordUF <- ordinate(ps.beta, "PCoA", "unifrac")
ordwUF <- ordinate(ps.beta, "PCoA", "wunifrac")
smpID <- sample_data(ps.beta)$MicroSam

# Keep first 2 vectors (latent variables, PCs) of each distance matrix
df <- rbind(data.frame(ordBC$vectors[,1:2], sample = smpID, method = 'BC'),
            data.frame(ordJC$vectors[,1:2], sample = smpID,method = 'Jaccard'),
            data.frame(ordUF$vectors[,1:2], sample = smpID,method = 'unifrac'),
            data.frame(ordwUF$vectors[,1:2], sample = smpID,method = 'wunifrac'))
# add sample_data info
df <- merge(df, data.frame(sample_data(ps.beta)), by.x= "sample",by.y= 'MicroSam')

# Calculate the mean for each level of MicroCatg
summary_df <- df %>%
  group_by(method,Microbiome,Category., MicroCatg) %>%
  summarise(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE)
  )
#Order Category. levels
summary_df$Category. <- factor(summary_df$Category., 
                               levels = c("Botanicals","Complex fiber","Gums","Oligosaccharrides",
                               "Waste Stream","b-glucan", "FBB16","FBB0"))
#Plot without arrows
ggplot(data = summary_df, aes(Axis.1, Axis.2, color = Microbiome, shape = Category.)) + 
  geom_point(size=3) + 
  facet_wrap(~method, scales = 'free') +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(2:7,1,19))  # Define shape values

#Set up FBB0 as the starting points
summary_df0<-summary_df
fbb0_df <- summary_df0 %>%
  filter(Category. == "FBB0") %>%
  select(method, Microbiome, Axis.1, Axis.2) %>%
  mutate(Axis.1_start = Axis.1, Axis.2_start = Axis.2) %>%
  select(-Axis.1, -Axis.2)  # Remove old columns to avoid duplicates

summary_df0 <- summary_df0 %>%
  left_join(fbb0_df, by = c("Microbiome", "method"))  # Match on method & Microbiome

#Plot with arrows
ggplot(summary_df0, aes(Axis.1, Axis.2, color = Microbiome, shape = Category..x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(2:7,1,19)) +  # Adjust shapes as needed
  labs(shape = "Category", color = "Microbiome") +  # Rename legends
  theme_minimal()
####################

#Subset data frame with only Weighted Unifrace Method
df0<-df[,c("Axis.1","Axis.2","Microbiome","MicroCatg","Category.","Sample","Highest.Dose.","method")]
df0<-df0[df0$method=="wunifrac",]

####Boxplot of Distances for Each Microbiome
#Set a Function to Plot the Distances to FBB0
####################
process_microbiome <- function(df0, Ki, mycolor) {
  library(dplyr)
  library(ggplot2)
  
  # Subset for FBB0
  df_FBB0 <- df0 %>% filter(Sample == "FBB0" & Microbiome == Ki)
  
  # Calculate the mean Axis.1 & Axis.2 values for FBB0
  FBB0_Axis1 <- mean(df_FBB0$Axis.1, na.rm = TRUE)
  FBB0_Axis2 <- mean(df_FBB0$Axis.2, na.rm = TRUE)
  
  # Remove FBB0 from the dataset
  df_Ki <- df0 %>% filter(Sample != "FBB0" & Microbiome == Ki)
  
  # Reset Axis.1 and Axis.2 based on the new coordinate axis origin (FBB0_Axis1, FBB0_Axis2)
  df_Ki <- df_Ki %>%
    mutate(Axis.1 = Axis.1 - FBB0_Axis1,
           Axis.2 = Axis.2 - FBB0_Axis2)
  
  # Calculate distances based on Axis.1 and Axis.2 values
  df_Ki <- df_Ki %>%
    mutate(distance = sqrt(Axis.1^2 + Axis.2^2))
  
  # Order Category levels
  df_Ki$Category. <- factor(df_Ki$Category., 
                            levels = c("Botanicals", "Complex fiber", "Gums", "Oligosaccharrides",
                                       "Waste Stream", "b-glucan", "FBB16"))
  
  # Box plot: Distance by Category
  p1 <- ggplot(df_Ki, aes(x = Category., y = distance, fill = Category.)) + 
    geom_boxplot(outlier.shape = 16,   
                 outlier.color = "black", 
                 outlier.size = 4) +
    labs(x = "Category", y = "Distance", 
         title = paste(Ki, "Category-Wise Distribution of Weighted UniFrac Distances to FBB0")) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1")
  
  print(p1)
  
  # Order Sample levels
  df_Ki$Sample <- factor(df_Ki$Sample, 
                            levels = c("beetroot powder", "BetaVia Pure WD",                     
                                       "carrot powder", "Cellulose",                           
                                       "dried apple pomace", "dried spinich powder",                
                                       "ginger powder", "Ground Beet Pulp","guar gum",                  
                                       "HiSmooth flax seed fiber", "Jerusalum articohoke powder (inulin)",
                                       "kappa carageenan ", "lemon powder",                       
                                       "locust bean gum", "Micronized cassava fiber",            
                                       "Mushroom chitosan", "Oat Bagasse hydrolysate-(B)",         
                                       "Oat Bagasse hydrolysate-(PBV)", "Oat Bagasse-Control",                 
                                       "Olive Flour", "pomegranate powder",                  
                                       "Psyllium Husk","xanthan gum","FBB16"))
  
  # Box plot: Distance by Sample
  p2 <- ggplot(df_Ki, aes(x = Sample, y = distance, fill = Sample)) + 
    geom_boxplot(outlier.shape = 16,   
                 outlier.color = "black", 
                 outlier.size = 1) +
    labs(x = "Sample", y = "Distance", 
         title = paste(Ki, "Distribution of Weighted UniFrac Distances from Samples to FBB0")) +
    theme(axis.text.x = element_text(angle = 30, hjust = 0.5, size = 7),
          legend.key.size = unit(17, "pt")) +
    scale_fill_manual(values = mycolor) +
    guides(fill = guide_legend(ncol = 1))
  
  print(p2)
  
  # Return the processed dataframe
  return(df_Ki)
}

####################

#Plot the Category and Sample Distances to FBB0
df_K3 <- process_microbiome(df0, "K3", mycolor)
df_K5 <- process_microbiome(df0, "K5", mycolor)
df_K7 <- process_microbiome(df0, "K7", mycolor)
df_K10 <- process_microbiome(df0, "K10", mycolor)


#Write a function to do Amova test to calculate the Fst and p value for overall and each pair of sample
#################
perform_amova_for_categories <- function(df_Ki) {
  
  # Compute overall AMOVA for all categories
  distance_matrix <- dist(df_Ki$distance)
  overall_amova <- adonis2(distance_matrix ~ Category., data = df_Ki, permutations = 999)
  
  overall_result <- data.frame(
    Category1 = "Overall",
    Category2 = "",
    Fs_value = overall_amova$F[1],
    p_value = overall_amova$`Pr(>F)`[1]
  )
  
  # Generate all unique pairs of categories
  category_pairs <- expand.grid(Category1 = unique(df_Ki$Category.), 
                                Category2 = unique(df_Ki$Category.)) %>%
    mutate(across(everything(), as.character)) %>%  # Convert factor to character
    filter(Category1 < Category2)  # Remove duplicate pairs (A,B) vs (B,A)
  
  # Initialize an empty results dataframe
  amova_results <- data.frame(Category1 = character(), Category2 = character(), 
                              Fs_value = numeric(), p_value = numeric())
  
  # Loop through each pair and perform AMOVA
  for (i in 1:nrow(category_pairs)) {
    c1 <- category_pairs$Category1[i]
    c2 <- category_pairs$Category2[i]
    
    # Subset only the two categories being compared
    df_pair <- df_Ki %>% filter(Category. %in% c(c1, c2))
    
    # Compute distance matrix for the pair
    distance_matrix_pair <- dist(df_pair$distance)
    
    # Run AMOVA for the pair
    amova_result <- adonis2(distance_matrix_pair ~ Category., data = df_pair, permutations = 999)
    
    # Extract Fs and p-value
    Fs_value <- amova_result$F[1]
    p_value <- amova_result$`Pr(>F)`[1]
    
    # Append results to dataframe
    amova_results <- rbind(amova_results, data.frame(Category1 = c1, 
                                                     Category2 = c2, 
                                                     Fs_value = Fs_value, 
                                                     p_value = p_value))
  }
  
  # Combine overall and pairwise results
  final_results <- rbind(overall_result, amova_results)
  
  return(final_results)
}

#################

#Calcualte Fst and p value for each pair of categories for each microbiome
df_K3_amova_categories <- perform_amova_for_categories(df_K3)
df_K5_amova_categories <- perform_amova_for_categories(df_K5)
df_K7_amova_categories <- perform_amova_for_categories(df_K7)
df_K10_amova_categories <- perform_amova_for_categories(df_K10)

write.csv(df_K3_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K3_amova_categories.csv")
write.csv(df_K5_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K5_amova_categories.csv")
write.csv(df_K7_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K7_amova_categories.csv")
write.csv(df_K10_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/haohaoK10_amova_categories.csv")

#Write a function to do Amova test to calculate the Fst and p value for overall and each pair of sample
#################
perform_amova_for_samples <- function(df_Ki) {
  
  # Compute overall AMOVA for all samples
  distance_matrix <- dist(df_Ki$distance)
  overall_amova <- adonis2(distance_matrix ~ Sample, data = df_Ki, permutations = 999)
  
  overall_result <- data.frame(
    Sample1 = "Overall",
    Sample2 = "",
    Fs_value = overall_amova$F[1],
    p_value = overall_amova$`Pr(>F)`[1]
  )
  
  # Generate all unique pairs of samples
  sample_pairs <- expand.grid(Sample1 = unique(df_Ki$Sample), 
                              Sample2 = unique(df_Ki$Sample)) %>%
    mutate(across(everything(), as.character)) %>%  # Convert factor to character
    filter(Sample1 < Sample2)  # Remove duplicate pairs (A,B) vs (B,A)
  
  # Initialize an empty results dataframe
  amova_results <- data.frame(Sample1 = character(), Sample2 = character(), 
                              Fs_value = numeric(), p_value = numeric())
  
  # Loop through each pair and perform AMOVA
  for (i in 1:nrow(sample_pairs)) {
    s1 <- sample_pairs$Sample1[i]
    s2 <- sample_pairs$Sample2[i]
    
    # Subset only the two samples being compared
    df_pair <- df_Ki %>% filter(Sample %in% c(s1, s2))
    
    # Compute distance matrix for the pair
    distance_matrix_pair <- dist(df_pair$distance)
    
    # Run AMOVA for the pair
    amova_result <- adonis2(distance_matrix_pair ~ Sample, data = df_pair, permutations = 999)
    
    # Extract Fs and p-value
    Fs_value <- amova_result$F[1]
    p_value <- amova_result$`Pr(>F)`[1]
    
    # Append results to dataframe
    amova_results <- rbind(amova_results, data.frame(Sample1 = s1, 
                                                     Sample2 = s2, 
                                                     Fs_value = Fs_value, 
                                                     p_value = p_value))
  }
  
  # Combine overall and pairwise results
  final_results <- rbind(overall_result, amova_results)
  
  return(final_results)
}

#################

#Calcualte Fst and p value for each pair of samples for each microbiome
df_K3_amova_sample <- perform_amova_for_samples(df_K3)
df_K5_amova_sample <- perform_amova_for_samples(df_K5)
df_K7_amova_sample <- perform_amova_for_samples(df_K7)
df_K10_amova_sample <- perform_amova_for_samples(df_K10)

write.csv(df_K3_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K3_amova_sample.csv")
write.csv(df_K5_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K5_amova_sample.csv")
write.csv(df_K7_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K7_amova_sample.csv")
write.csv(df_K10_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K10_amova_sample.csv")


#######Heatmap
####################
#glommate to Genus Level
ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])

#Plotting
heat.sample <- plot_taxa_heatmap(ps.genus,
                                 subset.top = 50,
                                 VariableA = "Category.",
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                 transformation = "log10")

#Calculate the average read counts for each category
################
# Extract OTU (abundance) table and metadata from phyloseq
otu_table_df <- as.data.frame(otu_table(ps.genus))
meta_df <- as.data.frame(sample_data(ps.genus))

# Add category information to the OTU table
otu_table_df$Category <- meta_df$Category.[match(rownames(otu_table_df), rownames(meta_df))]

# Aggregate by category (mean abundance per category)
otu_table_avg <- otu_table_df %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Convert back to matrix and ensure row names are correct
otu_table_avg <- column_to_rownames(otu_table_avg, var = "Category") %>%
  as.matrix()

# Convert back to a phyloseq OTU table
otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)

# Ensure tax_table matches the taxa in otu_table_avg_phylo
common_taxa <- intersect(rownames(tax_table(ps.genus)), colnames(otu_table_avg_phylo))
tax_table_filtered <- tax_table(ps.genus)[common_taxa, ]

# Create a new sample metadata table
meta_avg_phylo <- sample_data(data.frame(Category = rownames(otu_table_avg_phylo),
                                         row.names = rownames(otu_table_avg_phylo)))

# Set factor levels for Category in the desired order
meta_avg_phylo$Category <- factor(meta_avg_phylo$Category, 
                                  levels = c("Botanicals", "Complex fiber", "Gums",
                                             "Oligosaccharrides","Waste Stream",
                                             "b-glucan","FBB16","FBB0"))  # Replace with actual category names

# Create a new phyloseq object with matching OTU and tax tables
ps.genus.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)

# Generate heatmap of top 50 taxa, now grouped by category
heat.category <- plot_taxa_heatmap(ps.genus.avg,
                                   subset.top = 50,
                                   VariableA = "Category",
                                   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                   transformation = "log10")

# Plot the heatmap
heat.category

####################

#################


