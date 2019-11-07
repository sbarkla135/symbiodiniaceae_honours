####Sophie's code for processing bacterial 16S rRNA gene sequences from Symbiodiniaceae samples.
###Samples are subset into retentate and filtrate samples, as well as by Symbiodiniaceae strain.
###See subset samples at the end for full list of all subsets of samples created for this analysis.

##First, load the libraries you will need for analysis. Alternatively you can do this just before you need to use it
##Sometimes packages will not load simultaneously, and you will need to dismount some packages in order for others to run. If this
##happens, try restarting the program as sometimes this fixes the problem!

#################################################1. Load libraries.################################################
library(ape)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(RVAideMemoire)
library(indicspecies)

library(RColorBrewer)
library(gapminder)
library(microbiome)
library(BiocManager)

library(DESeq2)

################################################2. Import files from QIIME.######################################################

# The dataset is a phyloseq object!

#Importing the data in R as phyloseq object 

#Set working directory (where you have stored your QIIME output files)
setwd("~/Documents/Honours/R/final_R_files")

#Highlight the full set of each code to ensure it all runs together.

# Read in OTU table
otuData <- read.table(
  "ASV_unfiltered.csv",
  header = TRUE,
  sep = ",",
  row.names = 1
)


# Read in taxonomy table
taxData <- read.table(
  "taxonomy.csv",
  sep = ",",
  fill = TRUE,
  row.names = 1
)


# Read in metadata file
metaData <- import_qiime_sample_data('Metadata.txt')
#metadata file needs to be in a particular format, check on the google docs plugin called 'keemei' to ensure it is compatible with R.


# Add levels of taxonomy to taxonomy table
colnames(taxData) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


# Convert OTU and taxonomy tables from data frames to matrices for compatibility with phyloseq
otuDataMat <- as.matrix(otuData)
taxDataMat <- as.matrix(taxData)


# Import tree file as a phylosq object
tree <- read_tree("tree.nwk")


# Combine OTU, taxonomy, and metadata files into a phyloseq object
phy <- phyloseq(
  otu_table(
    otuDataMat, 
    taxa_are_rows = T
  ),
  tax_table(taxDataMat),
  sample_data(metaData)
)


# Merge the phyloseq objects into a single object
phy <- merge_phyloseq(phy, tree)


# Delete the matrices as they are not used in further analyses
rm(otuDataMat)
rm(taxDataMat)

#################################################3. Filtering################################################
## Filtering of unwanted taxa and samples 


# Remove chloroplast and mitochondria

phy <- subset_taxa(phy, Domain == "Bacteria" & Family != "mitochondria" & Class != "Chloroplast")

# Remove unneeded samples (Mocks)
phy <- subset_samples(phy, !(Sample_ID %in% c("Pos_control_Fwd19")))

# Function to remove rows and columns in
# a matrix or dataframe that sum to zero
# KillZeroRCs <- function(x) {
#   
#   x[ which( rowSums(x) != 0) , ] -> x
#   x[ , which( colSums(x) != 0) ] -> x
#   return(x)
#   }
#or

phy <- prune_taxa((taxa_sums(phy) > 0), phy)

#4. Decontamination of samples, removing contaminants identified by higher prevalence in control samples.
#You will need a column in your metadata file that identifies if a sample is a control or not for this code.
#in this code the column name is "Neg", and samples are listed as 'TRUE' or 'FALSE' in the metadata file.

library(decontam)

# identify contaminants
consList <- isContaminant(seqtab = phy, neg = "Neg", method = "prevalence")



# pull out the names of contaminants
cons <- rownames(consList)[consList$contaminant=="TRUE"]



# get info on the contaminants
# uncomment the following code to run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU
# then combine the consPer.csv and taxonomy.csv file data


# 
# subset the non neg control samples
vvv <- subset_samples(phy1, Neg == "FALSE")
# merge the samples
yyy <- merge_samples(vvv, "ID", fun = sum)
# transform counts to percentages
yyy <- transform_sample_counts(yyy, function(x) 100 * x/sum(x))
# extract the cons percentage data
zzz <- prune_taxa(x = yyy, taxa = cons)
# write otu table to dataframe
xxx <- data.frame(t(zzz@otu_table))
# write xxx to csv
write.csv(x = xxx, row.names = TRUE, file = "consPer.csv")
# subset the contaminant ASVs
probPhyCons <- prune_taxa(probPhy, taxa = cons)
probPhyCons <- prune_taxa(phy1, taxa = cons)
# write the contaminants to a file for reference
contaminants <- write_phyloseq(x = probPhyCons, type = 'TAXONOMY')

library(microbiome)

# remove the contaminants from the main phyloseq file
phy <- remove_taxa(phy, taxa = cons)

################################################# 4. Rarefaction ################################################
# Random subsampling of samples to ensure even depth is analysed from each sample. If asymptote is not reached for some samples, this
# means that the full diversity is not captured in your analysis.

library(vegan)

# # check counts per sample
sort(sample_sums(phy))
# # lowest value is e.g. 12345 for the samples I care about i.e. not the negative controls


# generate rarefied object
#sample.size is lowest value of counts per sample
#022 and 055 samples have really low reads, go to next highest one after blanks. 3722 is from sample 086AF
phyRare <- rarefy_even_depth(phy, sample.size = 10000, rngseed = 1)
# some low-count samples may be removed, but probably just the negative controls, so that's ok

##Console spit-out###

# `set.seed(1)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(1); .Random.seed` for the full vector
# ...
# 24 samples removedbecause they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 
#   
#   X22A_Fwd16X22AF_Fwd14X22B_Fwd16X22C_Fwd16X22CF_Fwd14
# ...
# 54OTUs were removed because they are no longer 
# present in any sample after random subsampling

#################


##check samples left after rarefying
# this allows you to check what samples were lost due to rarefaction.
sort(sample_sums(phyRare))


##### Generate rarefaction curve.
# "otuData" is the raw OTU table you imported from QIIME and use to create the a phyloseq object.
# Check rarecurve in the 'Help' tab to get info on 'step' and 'sample' parameters, but start with these values just to see what happens.

rarecurve(t(otuData), step = 100, sample = 10000, ylab = "Sample ASVs", xlab = "Sample reads", main = "Rarefaction Curve", cex = 0)


################################################# 5. Alpha Diversity ################################################

# Call phyloseq function
rich <- estimate_richness(phyRare2, split = TRUE, measures = NULL)

Sample_ID <- rownames(rich)

rich <- cbind(Sample_ID, rich) %>% as.data.frame()

write.table(rich, file = "Alphadiversity_allsamples.csv", sep = ",")



##t-test for p-values of different alpha diversity indices
ttest <- t(sapply(rich, function(x) unlist(t.test(x~sample_data(phyRare2)$Retentate_Filtrate)[c("estimate","p.value","statistic","conf.int")])))#to get alpha diversity stats
ttest



write.table(ttest, file = "Alphadiversity_ttest_withoutoutliers.csv", sep = ",") #output


# Adding the evenness measurement
rich$Evenness <- rich$Shannon / log(rich$Observed)

# Joining the metadata of the samples to the diversity indices

rich_fac <- left_join(sample_data(phyRare2), rich, by = "Sample_ID")

#Warning message:
#In class(x) <- c(subclass, tibble_class) :
#Setting class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an S4 object
#dont think this will affect anything but maybe check

#alpha diversity graph which shows all indexes: 

coral.sub4 <- subset_samples(phyRare2, SCF_RF %in% c("All"))


p <- plot_richness(phyRare2, "Retentate_Filtrate")

(p <- p + geom_boxplot(data = p$data, aes(x = Retentate_Filtrate, y = value, color = NULL), 
                       alpha = 0.1))

# Select samples of interest to visualise:
# do for all factors and all alpha-diversity indexes 
rich_fac_sub <- subset.data.frame(rich_fac, Retentate_Filtrate %in% c("Retentate", "Filtrate"))

ggplot(rich_fac_sub, aes(x = Retentate_Filtrate , y = Observed, fill = Retentate_Filtrate)) +
  geom_boxplot(position =position_dodge(width  =.8)) +
  theme_bw(base_size = 10)


######RUBY ALPHA DIVERSITY#########
p <- plot_richness(phyRare1, "Retentate_Filtrate")

#boxplots
(p <- p + geom_boxplot(data = p$data, aes(x = Retentate_Filtrate, y = value, color = NULL), 
                       alpha = 0.1))

(p <- p + geom_boxplot(data = p$data, aes(x = SCF_RF, y = Shannon, color = NULL), 
                       alpha = 0.1))

(p <- p + geom_boxplot(data = p$data, aes(x = Retentate_Filtrate, y = value, color = NULL), 
                       alpha = 0.1))


#bar graphs
(p <- p + geom_bar(data = p$data, aes(x = SymGenus, y = value, color = NULL), 
                   alpha = 0.1))

ggplot(rich_fac, aes(x = Retentate_Filtrate , y = Shannon, color = NULL)) +
  geom_bar(position =position_dodge(width  =.8)) +
  theme_bw(base_size = 10) + theme(axis.text.x=element_text(angle = 90))

ggplot(rich_fac, aes(x = Filter_Filtrate , y = Simpson, color = NULL)) +
  geom_bar(stat="identity") + 
  #geom_errorbar(aes)
  theme_bw(base_size = 10) + theme(axis.text.x=element_text(angle = 90)) +
  facet_wrap(vars(Filter_Filtrate), labeller = "label_value")

#if you want to reorder your samples in the bar graph 
#e.g. to have all retentate samples together, and all filtrate samples together
#change character to factor:
q$SymGenus_RF <- as.factor(q$SymGenus_RF)

#Reorder samples: write order after c() 
#write in order retentates next to each other and filtrates next to each other##
q$SymGenus_RF <- factor(q$SymGenus_RF, levels = c("A13_F", "B1_F", "C1a_F", "C1c_F","D1a_F", "D1b_F", 
                                                  "F1_F", "F5.1_F", "A13_R", "A3_R", "B1_R", "C1a_R", 
                                                  "C1b_R", "C1c_R", "D1a_R", "D1b_R", "F1_R"))


################################################# 6. Bar charts!!! ################################################
#Now you can finally visualise your data!

library(phyloseq)
library(ggplot2)
library(colorspace)
pal <- choose_palette()
set.seed(1)

# Load barchart colour pallete containing 500+ colours
source('~/Documents/Honours/R/Code/barColours1.R')
#This is just code that allows you to change the colours from the default rainbow settings which make it impossible to differentiate
#between objects. Can make your own if you go online and get different codes for colours in R.


#### Merge repeats/samples to have abundance go to 100 instead of 300###
#If your bar charts go over 100, you need to merge your samples so that each sample you want listed on the x-axis is merged together
#use subest_samples for this

# Transform counts to percentages (replace blanks with your phyloseq object)
phyObjectPer <- transform_sample_counts(blanks, function(x) 100 * x/sum(x))

# Collapse taxa to specified level of taxonomy (the example below is Family)
phyObjectPerFam <- tax_glom(phyObjectPer, "Genus")

# Melt phyloseq data
q <- psmelt(phyObjectPerFam)

#change character to factor:
#q$SymGenus_RF <- as.factor(q$SymGenus_RF)

#Reorder samples: write order after c() 
#write in order retentates next to each other and filtrates next to each other##
#q$SymGenus_RF <- factor(q$SymGenus_RF, levels = c("A13_F", "B1_F", "C1a_F", "C1c_F","D1a_F", "D1b_F", 
#                                                 "F1_F", "F5.1_F", "A13_R", "A3_R", "B1_R", "C1a_R", 
#                                                "C1b_R", "C1c_R", "D1a_R", "D1b_R", "F1_R"))

#"A13_F"  "A13_R"  "A3_R"   "B1_F"   "B1_R"   "C1a_F"  "C1a_R"  "C1b_R"  "C1c_F"  "C1c_R"  "D1a_F" 
#[12] "D1a_R"  "D1b_F"  "D1b_R"  "F1_F"   "F1_R"   "F5.1_F"

#lists all names in this level to get printout of names for reordering (above)
#levels(q$SymGenus_RF)

#facet_wrap to split a ggplot into two separate bar graphs in one graph - e.g. group retenetate and filtrates.
# use -> scales = "free_x" - in the facet_wrap line of code, to keep sample names listed only in their groups, and not twice (e.g. in both retentate and filtrate sections of graph)


# Generate ggplot2 object
p <- ggplot(q, aes_string(x = "ID", y = "Abundance", fill = "Genus"))

# Output ggplot2 object with modifications

p +  
  geom_bar(stat = 'identity', position = "stack") +
  ylab("Relative Abundance") +
  xlab("ID")+
  theme_bw() +
  theme(legend.position = "right") +
  scale_fill_manual(values = barColours1) + # If barColours1 looks bad, try barColours2 or 3
  guides(fill = guide_legend(ncol = 4)) + 
  # Adjust the number of columns to suit
  theme(axis.text.x = element_text(angle = 90)) +
  #facet_wrap(vars(Retentate_Filtrate), labeller = "label_value", scales = "free_x")
  
  
  # Note: Some bars may not extend to 100 because many ASVs cannot be identified to the
  # taxonomic level you are analysing.
  
  # IMPORTANT !!! Your barchart is likely to be very big. It will probably be be too big to
  # display in the R Studio plots window. To preview and/or output, use the Export dropdown
  # menu in the plots window. Select PDF. Adjust the size by entering some large dimensions
  # e.g. 33 inches x 33 inches. Hit the preview button to preview. Hit the save/export button
  # to save/export.

################################################# 7. Beta Diversity #######################################################

# Good idea to set the seed before doing these kinds of analysis in order to generate the same outcome every time
set.seed(1)

# Calculating dissimilarity matrix:
DistBC.coral = phyloseq::distance(phyRare2, method = "bray") # Relative abundance; could also be Jaccard for presence/absence
dat.coral = as(sample_data(phyRare2), "data.frame")

# Performing the permutational analysis of variance on all groups:
# I am considering here the variable "Stage", which represents the life stage of coral

#for permanova on alpha diversity, read in alpha diversity table generated in alpha diversity code
alpha <- read.table(
  "alpha.metadata.csv",
  header = TRUE,
  sep = ",",
  row.names = 1
)

alpha.aov <- aov(Shannon~SymGenus*Retentate_Filtrate, data = alpha)
plot(alpha.aov)
summary(alpha.aov)


adonis2(DistBC.coral ~ SymGenus, strata = dat.coral$SymGenus, data = dat.coral)

adonis2(DistBC.coral ~ SymGenus + Retentate_Filtrate + SymGenus*Retentate_Filtrate, strata = dat.coral$SymGenus + Retentate_Filtrate + SymGenus*Retentate_Filtrate, data = dat.coral)

interaction.plot(alpha$SymGenus, alpha$Retentate_Filtrate, alpha$Observed)

###try nvabund### because permanova results for within strain between retentate and filtrate isnt significant?? dave warton
#permanova thros
####################################Sophie results####################################
######For strain

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = DistBC.coral ~ Strain, data = dat.coral, strata = dat.coral$Strain)
# Df SumOfSqs     R2      F Pr(>F)    
# Strain   12   7.7191 0.3707 2.3072  0.001 ***
#   Residual 47  13.1036 0.6293                  
# Total    59  20.8227 1.0000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

####For A13 filter/filtrate comparison
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 119
# 
# adonis2(formula = DistBC.coral ~ Filter_Filtrate, data = dat.coral, strata = dat.coral$Filter_Filtrate)
# Df SumOfSqs      R2     F Pr(>F)
# Filter_Filtrate  1  1.06372 0.88558 23.22    0.1
# Residual         3  0.13743 0.11442             
# Total            4  1.20115 1.00000  

######For B1 filter/filtrate
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 119
# 
# adonis2(formula = DistBC.coral ~ Filter_Filtrate, data = dat.coral, strata = dat.coral$Filter_Filtrate)
# Df SumOfSqs      R2      F Pr(>F)
# Filter_Filtrate  1  0.86672 0.66984 6.0866    0.1
# Residual         3  0.42720 0.33016              
# Total            4  1.29391 1.00000   

####C1a filter/filtrate
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 119
# 
# adonis2(formula = DistBC.coral ~ Filter_Filtrate, data = dat.coral, strata = dat.coral$Filter_Filtrate)
# Df SumOfSqs      R2      F Pr(>F)
# Filter_Filtrate  1  1.15538 0.93938 46.487    0.1
# Residual         3  0.07456 0.06062              
# Total            4  1.22994 1.00000  

##C1b not enough samples left after rarefying to perform this test

#########C1c filter/filtrate
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 719
# 
# adonis2(formula = DistBC.coral ~ Filter_Filtrate, data = dat.coral, strata = dat.coral$Filter_Filtrate)
# Df SumOfSqs     R2      F Pr(>F)
# Filter_Filtrate  1  0.74142 0.4981 3.9697    0.1
# Residual         4  0.74708 0.5019              
# Total            5  1.48850 1.0000  

#########D1a filter/filtrate
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 719
# 
# adonis2(formula = DistBC.coral ~ Filter_Filtrate, data = dat.coral, strata = dat.coral$Filter_Filtrate)
# Df SumOfSqs      R2      F Pr(>F)
# Filter_Filtrate  1  1.38536 0.85854 24.276    0.1
# Residual         4  0.22827 0.14146              
# Total            5  1.61364 1.00000   

#########D1b filter/filtrate
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 119
# 
# adonis2(formula = DistBC.coral ~ Filter_Filtrate, data = dat.coral, strata = dat.coral$Filter_Filtrate)
# Df SumOfSqs     R2      F Pr(>F)
# Filter_Filtrate  1 0.239310 0.9591 70.346    0.1
# Residual         3 0.010206 0.0409              
# Total            4 0.249516 1.0000  

#########F1 filter/filtrate
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 719
# 
# adonis2(formula = DistBC.coral ~ Filter_Filtrate, data = dat.coral, strata = dat.coral$Filter_Filtrate)
# Df SumOfSqs      R2      F Pr(>F)
# Filter_Filtrate  1  1.36849 0.86408 25.429    0.1
# Residual         4  0.21526 0.13592              
# Total            5  1.58375 1.00000  

###not enough samples left for F5.1 analysis after rarefying

#########For filtrate samples between strains
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = DistBC.coral ~ SymGenus, data = dat.coral, strata = dat.coral$SymGenus)
# Df SumOfSqs      R2      F Pr(>F)    
# SymGenus  7   7.4363 0.94826 41.889  0.001 ***
#   Residual 16   0.4058 0.05174                  
# Total    23   7.8420 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#######For retentate samples between strains
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = DistBC.coral ~ SymGenus, data = dat.coral, strata = dat.coral$SymGenus)
# Df SumOfSqs      R2      F Pr(>F)   
# SymGenus  8   3.0694 0.64742 2.7543  0.001 **
#   Residual 12   1.6716 0.35258                 
# Total    20   4.7409 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###this stat significance could be due to the outliers####

##retenatesamples2 - retentate with outliers removed
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = DistBC.coral ~ SymGenus, data = dat.coral, strata = dat.coral$SymGenus)
# Df SumOfSqs      R2      F Pr(>F)
# SymGenus  6  0.61951 0.37252 0.9894  0.533
# Residual 10  1.04353 0.62748              
# Total    16  1.66305 1.00000  


##two factor crossed (two fixed factors) model - ANOVA - page 231 quinn & keough 
##with outliers
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = DistBC.coral ~ SymGenus + Retentate_Filtrate + SymGenus * Retentate_Filtrate, data = dat.coral, strata = dat.coral$SymGenus + Retentate_Filtrate + SymGenus * Retentate_Filtrate)
# Df SumOfSqs      R2       F Pr(>F)    
# SymGenus                     9   7.7483 0.46540 11.5908  0.001 ***
#   Retentate_Filtrate           1   3.5562 0.21360 47.8778  0.001 ***
#   SymGenus:Retentate_Filtrate  6   3.2646 0.19608  7.3253  0.001 ***
#   Residual                    28   2.0797 0.12492                   
# Total                       44  16.6488 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##without outliers
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = DistBC.coral ~ SymGenus + Retentate_Filtrate + SymGenus * Retentate_Filtrate, data = dat.coral, strata = dat.coral$SymGenus + Retentate_Filtrate + SymGenus * Retentate_Filtrate)
# Df SumOfSqs      R2       F Pr(>F)    
# SymGenus                     8   6.5495 0.43756 14.6570  0.001 ***
#   Retentate_Filtrate           1   4.2636 0.28485 76.3327  0.001 ***
#   SymGenus:Retentate_Filtrate  5   2.7029 0.18058  9.6781  0.001 ***
#   Residual                    26   1.4523 0.09702                   
# Total                       40  14.9683 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

####two factor anova for alpha diversity
##observed alpha diversity
# > summary(alpha.aov)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# SymGenus                     8   3891   486.4   6.105 0.000188 ***
#   Retentate_Filtrate           1      1     1.5   0.019 0.892570    
# SymGenus:Retentate_Filtrate  5    979   195.8   2.457 0.059650 .  
# Residuals                   26   2072    79.7                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###simpson alpha diversity
# Df  Sum Sq Mean Sq F value   Pr(>F)    
# SymGenus                     8 0.09037 0.01130   3.331  0.00934 ** 
#   Retentate_Filtrate           1 0.20572 0.20572  60.660 2.91e-08 ***
#   SymGenus:Retentate_Filtrate  5 0.03794 0.00759   2.238  0.08066 .  
# Residuals                   26 0.08817 0.00339                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##shannon alpha diversity
#                            Df Sum Sq Mean Sq F value   Pr(>F)    
# SymGenus                     8  2.991   0.374   7.618 3.24e-05 ***
#   Retentate_Filtrate           1 13.127  13.127 267.438 3.36e-15 ***
#   SymGenus:Retentate_Filtrate  5  0.870   0.174   3.547   0.0141 *  
#   Residuals                   26  1.276   0.049                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

write.csv(as.data.frame(adonis2), file="adonis2-results.csv")

# Checking for homogeneity of variance:
# If the p-value is significiant (<0.05) -> variances are not homogeneous.
# However, this is still acceptable if our design is balanced
beta <- betadisper(DistBC.coral, sample_data(phyRare2)$Retentate_Filtrate)
disper.test = permutest(beta)
disper.test # Ok if p > 0.05

#p-value for strain type is 0.003 **
#for extraction.date p = 0.001***
#for experiement.stage p = 0.001***
#for repllicate.number p = 0.001***
#for temperature p = 0.001***

##Sophie results##

#FOR STRAIN#
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups    12 0.6969 0.058073 0.5742    999  0.856
# Residuals 47 4.7535 0.101138  

###For strain on all samples###
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
########## Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     8 0.8884 0.11105 0.8519    999  0.563
# Residuals 32 4.1712 0.13035     

######For association(retentate/filtrate) on all samples####
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.62372 0.62372 57.348    999  0.001 ***
#   Residuals 39 0.42417 0.01088                         
---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  ##From these results, on level of strain - there was homogenous variance, but not between retentate
  #and filtrate samples. Therefore pairwise should only be done on strain level
  
  # If homogeneity of variances is respected, we can do pairwise comparisons
  
  ##p-value for strain type in filtrate samples = 0.038
  
  testingBC.coral = pairwise.perm.manova(DistBC.coral, sample_data(phyRare1)$Filter_Filtrate,
                                         nperm=999, p.method = "BH")

testingBC.coral$p.value

#replicate          A        B        C        D         E  F
#B 1.000000       NA       NA       NA        NA NA
#C 1.000000 1.000000       NA       NA        NA NA
#D 1.000000 1.000000 1.000000       NA        NA NA
#E 1.000000 1.000000 1.000000 1.000000        NA NA
#F 0.112875 0.112875 0.112875 0.132300 0.1306667 NA
#N 0.112875 0.112875 0.112875 0.112875 0.1128750  1

##SOPHIE RESULTS##
#strain##
# > testingBC.coral$p.value
# B1     Blank PCR_blank       S22       S26       S49       S55       S65       S82
# Blank     0.4226182        NA        NA        NA        NA        NA        NA        NA        NA
# PCR_blank 1.0000000 1.0000000        NA        NA        NA        NA        NA        NA        NA
# S22       0.6830909 0.9887324 1.0000000        NA        NA        NA        NA        NA        NA
# S26       0.2705625 0.4862000 0.8396471 0.4217778        NA        NA        NA        NA        NA
# S49       0.2860000 0.4747241 0.9653913 0.3533878 0.2775882        NA        NA        NA        NA
# S55       0.3493750 0.9887324 1.0000000 0.6000000 0.3319149 0.3029302        NA        NA        NA
# S65       0.2705625 0.2705625 0.5162857 0.2775882 0.2705625 0.3310667 0.2860000        NA        NA
# S82       0.2705625 0.4638214 1.0000000 0.3319149 0.2705625 0.4747241 0.2860000 0.4118400        NA
# S86       0.2652000 0.2705625 0.4217778 0.2705625 0.2652000 0.2705625 0.2705625 0.2705625 0.2705625
# S89       0.2705625 0.4798983 1.0000000 0.4217778 0.2705625 0.2705625 0.2860000 0.2705625 0.2705625
# S92       0.2705625 0.3027143 0.5277188 0.3173182 0.2705625 0.2652000 0.2705625 0.2705625 0.2705625
# Water     0.4217778 0.8149254 1.0000000 0.5032258 0.2886000 0.3027143 0.5032258 0.2705625 0.2860000
# S86       S89       S92
# Blank            NA        NA        NA
# PCR_blank        NA        NA        NA
# S22              NA        NA        NA
# S26              NA        NA        NA
# S49              NA        NA        NA
# S55              NA        NA        NA
# S65              NA        NA        NA
# S82              NA        NA        NA
# S86              NA        NA        NA
# S89       0.2652000        NA        NA
# S92       0.2652000 0.2705625        NA
# Water     0.2705625 0.2705625 0.2705625


print(testingBC.coral)

#Pairwise comparisons using permutation MANOVAs on a distance matrix 

#filter_filtrate

#Filter
#Filtrate 0.001 
#999 permutations

#csv output

write.table(testingBC.coral[["p.value"]], file="permanova.csv", sep=",")

# PCoA plot (you can do nMDS instead)
coral.sub <- subset_samples(phyRare1, Retentate_Filtrate %in% c("Retentate", "Filtrate"))
coral.sub2 <- subset_samples(retentatesamples, Retentate_Filtrate %in% c("Retentate"))
coral.sub3 <- subset_samples(filtratesamples, Retentate_Filtrate %in% c("Filtrate"))
coral.sub4 <- subset_samples(class4.clean, experiment.stage %in% c("before", "after"))

#do this with wunifrac as well
class.pcoa <- ordinate(
  physeq = phyRare2, 
  method = "PCoA", 
  distance = "bray")

class2.pcoa <- ordinate(
  physeq = retentatesamples, 
  method = "PCoA", 
  distance = "bray")

class3.pcoa <- ordinate(
  physeq = retentatesamples2, 
  method = "PCoA", 
  distance = "bray")

class4.pcoa <- ordinate(
  physeq = filtratesamples, 
  method = "PCoA", 
  distance = "bray")

plot_ordination(
  physeq = phyRare2,
  ordination = class.pcoa,
  color = "SymGenus",
  shape = "Retentate_Filtrate",
  axes = 1:3,
  title = "PCoA of bacterial communities in all samples from Bray-Curtis") +
  geom_point(aes(color = SymGenus), size = 8) +
  theme_bw(base_size = 20)

plot_ordination(
  physeq = retentatesamples ,
  ordination = class2.pcoa,
  color = "Sample_ID",
  label = "Sample_ID",
  shape = "Filter_Filtrate",
  axes = 1:3,
  title = "PCoA of bacterial communities in retentate samples from Bray-Curtis") +
  geom_point(aes(color = Sample_ID), size = 8) +
  theme_bw(base_size = 20)

plot_ordination(
  physeq = phyRare,
  ordination = class3.pcoa,
  color = "SymGenus",
  shape = "Filter_Filtrate",
  axes = 1:3,
  title = "PCoA of bacterial communities in select retentate samples from Bray-Curtis") +
  geom_point(aes(color = SymGenus), size = 8) +
  theme_bw(base_size = 20)

plot_ordination(
  physeq = filtratesamples ,
  ordination = class4.pcoa,
  color = "SymGenus",
  shape = "Retentate_Filtrate",
  axes = 1:3,
  title = "PCoA of bacterial communities in filtrate samples") +
  geom_point(aes(color = SymGenus), size = 8) +

  
################################################# 8. DESEQ #######################################################

##this is done on unrarefied dataset - original phyloseq object you created
#You will need to have no zeros in your dataset - which you likely do. So you can just add 1 to each read, which barely effects data
#and makes the code work :D - do this as below (make a new phyloseq with the new Otu file)
#make sure you're naming the new OTU file and phyloseq to a new name, so you keep the original as well

#adding 1 to each count
otuData1 <- otuData+1

#now remerge phyloseq object- do not reimport the otu table, but reform into matrix then
#merge phyloseq object

# Convert OTU and taxonomy tables from data frames to matrices for compatibility with phyloseq
otuDataMat1 <- as.matrix(otuData1)
taxDataMat <- as.matrix(taxData)


# Combine OTU, taxonomy, and metadata files into a phyloseq object
phy1 <- phyloseq(
  otu_table(
    otuDataMat1, 
    taxa_are_rows = T
  ),
  tax_table(taxDataMat),
  sample_data(metaData)
)


# Merge the phyloseq objects into a single object
phy1 <- merge_phyloseq(phy1, tree)


# Delete the matrices as they are not used in further analyses
rm(otuDataMat1)
rm(taxDataMat)

#now let's get on with DeSeq...
# I then removed the samples I wanted to with subset_samples (e.g. outliers and controls), which is why I am using 'phy3'

#genuslevelphy <-  tax_glom(phy1, "Genus")

coral.subset <- subset_samples(phy3, Retentate_Filtrate %in% c("Retentate", "Filtrate"))


data.deseq2 <- phyloseq_to_deseq2(coral.subset, ~ Retentate_Filtrate) # converting in good format



dds <- DESeq(data.deseq2, test="Wald", fitType="parametric")


colData(dds)


# alpha <- 0.01 # setting significance threshold value
# sigtab <- data.deseq2.res[which(data.deseq2.res$padj < alpha), ] # putting in a table significant taxa
# sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(otuData1)[rownames(sigtab), ], "matrix"))


res <- results(dds, tidy=TRUE, contrast=c("Retentate_Filtrate", "Retentate", "Filtrate")) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
res

#p-adj is better than p-value when you have a lot of calculations, as the statistical significance can get skewed. 

## 1:15 means you are using the first 15 rows, so the top 15 most statistically significantly different. You can change this.
## you can check the file 'res' to check the p-values of each one, to help decide how many you want to display.
goi <- res$row[1:15]
stopifnot(all(goi %in% names(dds)))
goi



write.csv(as.data.frame(res), file="deseq-top15-PVALUE-nocontaminants.csv")



tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(Species, Abundance, (ncol(.)-length(goi)+1):ncol(.))



tcounts %>% 
  select(Row.names, Retentate_Filtrate, SymGenus_RF, Species, Abundance) %>% 
  head %>% 
  knitr::kable()



ggplot(tcounts, aes(Retentate_Filtrate, Abundance, fill=Retentate_Filtrate)) + 
  geom_boxplot() + 
  facet_wrap(~Species, scales="free") + 
  labs(x="Retentate_Filtrate", 
       y="Abundance (log normalized counts)", 
       fill="(group)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

q$SymGenus_RF <- as.factor(q$SymGenus_RF)
q$SymGenus_RF <- factor(q$SymGenus_RF, levels = c("A13_F", "B1_F", "C1a_F", "C1c_F","D1a_F", "D1b_F", 
                                                  "F1_F", "F5.1_F", "A13_R", "A3_R", "B1_R", "C1a_R", 
                                                  "C1b_R", "C1c_R", "D1a_R", "D1b_R", "F1_R"))

ggplot(tcounts, aes(x = ID, y = Species, size = Abundance, color = Retentate_Filtrate)) +
  geom_point() +
  scale_size_area(max_size = 6, breaks = c(1, 5, 10, 20, 40, 70, 90)) + # determines size of bubbles
  #coord_fixed(ratio = 0.9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(vars(Retentate_Filtrate), labeller = "label_value", scales = "free_x")

  
#Fin. Email me @ sophiebarkla@gmail.com for any queries :)