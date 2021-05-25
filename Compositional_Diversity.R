install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
remotes::install_github("jbisanz/qiime2R")
BiocManager::install("phyloseq")
install.packages("tidyverse")
install.packages("ggplot2")
BiocManager::install("microbiome")
install.packages("vegan")
# Loading feature tables as subsettable objects 
library(qiime2R)
library(phyloseq)
#Dataframe management, plotting
library(tidyverse)
library(ggplot2)
#Transformation and summary function
library(microbiome); packageVersion("microbiome")
#Species richness
library(breakaway)
# Statistical testing for diversity diifferences*
library(vegan)
setwd("~/Documents/Qiime2-NFLD_AZMP/18S/Trim-350")

# If you used QIIME2 you have qiime artefacts (.qza) for feature tables (ASVs x Samples), a phylogenetic tree, and metadata.
# qiime2R allows you to load all this data together, as a subsettable phyloseq object (S4 class)
physeq<-qza_to_phyloseq(
  features="deblur_output/deblur_table_final.qza",
  tree="tree_out/rep_seqs_final_aligned_masked_tree_rooted.qza",
  "taxa/classification.qza",
  metadata = "Metadata.txt"
)

#Check out the data
physeq %>% sample_data

#The independent variable used in this dataset is "Transect", you will have your own

#Summary function for the ASV data (check out ?summarize_phyloseq)
summarize_phyloseq(physeq)
# Singletons are ASVs observed once. The frequency counts of an ASV table are an informative property.
#It can be used to estimate unobserved diversity and apply some classic statistical
# methods for incomplete biological samples, analogous to accounting for unequal sampling efforts. 

Diatoms <- physeq %>% 
  subset_taxa(Phylum == "Diatomea")
otu_table(Diatoms)
tax_table(Diatoms)
sample_data(Diatoms)

############################### This section works on raw datasets: No transformations, just ASV counts ############################### 

#The alpha-diversity analysis is built through Amy Willis : https://github.com/adw96/stamps2018/blob/master/estimation/diversity-lab.R

#Amy has several papers directly relevant to diversity analysis in microbiome data
#Breakaway paper: Willis, A., & Bunge, J. (2015). Estimating diversity via frequency ratios: Estimating Diversity via Ratios. Biometrics, 71(4), 1042–1049. https://doi.org/10.1111/biom.12332.
#Willis, A., Bunge, J., & Whitman, T. (2015). Inference for changes in biodiversity. ArXiv:1506.05710 [q-Bio, Stat]. http://arxiv.org/abs/1506.05710.
#Willis, A., Bunge, J., & Whitman, T. (2017). Improved detection of changes in species richness in high diversity microbial communities. Journal of the Royal Statistical Society: Series C (Applied Statistics), 66(5), 963–977. https://doi.org/10.1111/rssc.12206.
#Willis, A. D. (2019). Rarefaction, Alpha Diversity, and Statistics. Frontiers in Microbiology, 10, 2407. https://doi.org/10.3389/fmicb.2019.02407.

#Unequal sequencing depths are a challenge for comparing diversity within (alpha) and between (beta) samples.
#breakaway is a package to estimate the number of missing taxa based on the frequency ratios of ASVs (i.e., the number of taxa observed once =singleton, observed twice = doubleton, etc.)
#The breakaway model estimates unobserved Genera (in this case) by fitting a series of non-linear regression models (e.g., Poisson, Negative Binomial), choosing a best fit
#to consecutive taxa frequency ratios (singletons, doubletons, tripletons, etc.) and predicting
#the number of taxa unobserved (i.e., The number of groups with a frequency of zero). 
#A more intuitive explanation is that if relatively few taxa were rarely observed, 
#then most of the diversity was probably sampled (More often the case with high sequencing depths). 
#thus standard errors for unobserved taxa will be small  and the observed richness is reliable, 
# and vice-versa. So, richness estimates are "adjusted" for sequencing depth by showing high SE 
# where there is high singleton frequency, more often the case in shallow depth samples. Phew.  

#Microbiome data does contain spurious singletons, and there is an alternative to breakaway (breakaway_nof1) 
# that disregards singletons when modelling frequency ratios, shown further down.

#Let's begin with the dataset in the phyloseq object
#If you want to agglomerate the taxonomy to a genus level, e.g., if you're interested in genus-level
#richness instead of species-level (especially because many ASVs are unclassified at a species level)
Genus <- physeq %>% 
  tax_glom("Genus")

#Inspect the observed richness in the dataset (The)
observed <- sample_richness(Genus)
summary(observed)
#Plot observed richness by some variable in sample_data you viewed above
plot(observed, Genus, color="Transect")

#Run breakaway model on Genus level object
estimated_richness <- breakaway(Genus)
plot(estimated_richness, Genus, color="Transect")
#Now the visualization includes SE for richness based on ASV frequency ratios
summary(estimated_richness)
# The breakaway model includes richness estimate, SE, and the chosen models to the observed frequency ratios (Potentially different between samples)

#Next-gen. seqencing datasets can contain false singletons due to (potentially systematic) sequencing errors,
# This will inflate the number of rare taxa and skew richness estimates based on frequency ratios.
#Exploratory tool is to fit same breakaway models to richness ignoring the singletons.
estimated_richness_nof1 <- breakaway_nof1(Genus)
plot(estimated_richness_nof1, Genus, color="Transect")
#SE are likely to be higher in some samples. Interesting. 

#Important to know that by reducing your dataset with singleton removal, Amy warns the models are empirically less "stable".
#Revisit the summarise_phyloseq above and consider singleton frequency (i.e., sparisty)

#Although I am unconvinced that hypothesis testing for diversity differences in environmental microbiome data has much weight,
# where most of the methods are from human clinical samples with much, much larger effect sizes, 
# it can be done.
bt <- betta(summary(estimated_richness)$estimate,
            summary(estimated_richness)$error,
            make_design_matrix(Genus, "Transect"))
view(bt$table)

# If you're interested to take this further, checkout the extension of breakaway to evenness (DivNet). 
#Here is nice starting point (https://adw96.github.io/breakaway/articles/intro-diversity-estimation.html)
# This will help you inspect the model fit for breakaway and test your data further

############################### This section works on transformed datasets: clr transform ############################### 

#This section touches beta-diversity: Calculating a compositionally-valid distance metric between 
#centered-log-ratio transformed ASVs. The "Aitchinson's" distance is that metric here, which is essentially
# a true euclidean distance (As opposed to dissimilarity e.g., Bray-curtis) of the clr data which is 
# sub-compositionally coherent (i.e., it respects the relative nature of the data).
# The beta-diversity can be represented in high-dimensional space to look for patterns- ordination can reveal continuous patterns in high dimensions.
# Because the clr-transform allows a euclidean distance, you can apply a Principal component analysis (PCA)
# as an ordination technique.

#clr-transformation on ASVs
physeq <- microbiome::transform(physeq, "clr")  

#Ordinate via redundancy analysis (RDA) which can be analyzed as a PCA
ord_clr <- phyloseq::ordinate(physeq, "RDA", distance = "euclidean")

#These are some plots I've used
#Scree plot (how much variation is explained in high dimensions by each axis)
scree <- phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "red") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")+
  theme_bw(28)+
  theme(axis.text.x.bottom = element_text(angle = 45, hjust=1))
scree

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
PCA <- phyloseq::plot_ordination(physeq, ord_clr, color="Transect") + 
  geom_point(size = 5) +
  coord_fixed(clr2 / clr1) +
  theme_bw(28)+
  stat_ellipse()+
  theme(legend.position = "bottom")
PCA

#Statistically test whether there is more variation between groups than within them (groups are transects in this case)

#Extract the compositional distances (i.e., measure of beta-diversity)
dis_clr <- phyloseq::distance(physeq, method  = "euclidean")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq))

# Here, the clr-distances is a function of transect variable
perm <- vegan::adonis(dis_clr ~ sampledf$Transect, data = sampledf)
perm

#There are lots methods you can use for alpha and beta diversity. The breakaway model is a good start
# for adjusting richness estimates by sequencing depths. Although very, very shallow samples seem likely to be incomparable
# As far as I understand, alpha diversity does not necessarily need to be compositional, because you are dealing with frequencies of presence/absence, not
# proportional changes. Beta-diversity should require compositional analysis, and a Aitchinson's distance
#combined with ordination and PERMANOVA are widely used.

#See 
#Calle, M. L. (2019). Statistical Analysis of Metagenomics Data. Genomics & Informatics, 17(1), e6. https://doi.org/10.5808/GI.2019.17.1.e6.
# for an overview