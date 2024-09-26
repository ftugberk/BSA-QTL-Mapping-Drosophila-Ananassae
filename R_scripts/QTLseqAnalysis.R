# Load necessary libraries
library(vcfR)
library(QTLseqr)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(clusterProfiler)
library(enrichR)
library(openxlsx)

# Set working directory
setwd('C:/Users/ftugb/OneDrive/Desktop/Master Thesis/R_analysis')

# Import data for tolerant and sensitive offspring and parents
poolTO <- read_tsv("fast_O.table")
poolSO <- read_tsv("slow_O.table")
poolTP <- read_tsv("fast_P.table")
poolSP <- read_tsv("slow_P.table")

# Rename columns for clarity
rename_columns <- function(df, prefix) {
  names(df)[names(df) == "sample01.AD"] <- paste0(prefix, ".AD")
  names(df)[names(df) == "sample01.DP"] <- paste0(prefix, ".DP")
  names(df)[names(df) == "sample01.GT"] <- paste0(prefix, ".GT")
  return(df)
}

poolTO <- rename_columns(poolTO, "TolerantOffspring")
poolSO <- rename_columns(poolSO, "SensitiveOffspring")
poolTP <- rename_columns(poolTP, "TolerantParent")
poolSP <- rename_columns(poolSP, "SensitiveParent")


offspring_pools <- full_join(poolSO, poolTO)
write_tsv(offspring_pools, "offspring_pools.table")

# Define files, chromosome list, high and low bulk 
file <- "offspring_pools.table"
Chroms <- c("CM029943.2", "CM029944.2", "CM033063.1", "CM033064.1", "CM029945.2", "CM029946.2")
HighBulk <- "TolerantOffspring"
LowBulk <- "SensitiveOffspring"

# Import data from GATK
df <- importFromGATK(file = file, highBulk = HighBulk, lowBulk = LowBulk, chromList = Chroms)

# Visualize Depth and Reference Frequency Distributions prior to filtration

# Histogram of total depth
ptotaldepth1 <- ggplot(df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  theme_minimal() +
  xlim(0,1000)

# Histogram of reference allele frequency
pallelefreq1 <- ggplot(df) +
  theme_minimal() +
  geom_histogram(aes(x = REF_FRQ))


# Filter SNPs and observe depth and reference frequency distribution.

df_filt <- filterSNPs(SNPset = df, refAlleleFreq = 0.20, depthDifference = 100, maxTotalDepth = 400, verbose = TRUE)
df_nona <- na.omit(df_filt)

# Visualize data after filtering
ptotaldepth2 <- ggplot(df_nona) +
  theme_minimal() +
  geom_histogram(aes(x = REF_FRQ))

pallelefreq2 <- ggplot(df_nona) +
  theme_minimal() +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  xlim(0,1000)

pSNPindexhigh <- ggplot(df_nona) +
  theme_minimal() +
  geom_histogram(aes(x = SNPindex.HIGH))

pSNPindexlow <- ggplot(df_nona) +
  theme_minimal() +
  geom_histogram(aes(x = SNPindex.LOW))


### QTLseq and G prime analysis

# Run QTLseq analysis
qtl_off <- runQTLseqAnalysis(df_nona, windowSize = 1e6, popStruc = "RIL", bulkSize = 100, replications = 1000000, intervals = c(95,99))

#Change GenBank chromosome names into their actual names
qtl_off$CHROM <- factor(qtl_off$CHROM, levels = c("CM029943.2", "CM029944.2", "CM033063.1", "CM033064.1", "CM029945.2", "CM029946.2"), labels = c("Chr 2L", "Chr 2R", "Chr 3L", "Chr 3R", "Chr XL", "Chr XR"))

# Run G' analysis
gprime <- runGprimeAnalysis(df_nona, windowSize = 1e6, outlierFilter = "deltaSNP", filterThreshold = 0.05)
#Change GenBank chromosome names into their actual names
gprime$CHROM <- factor(gprime$CHROM, levels = c("CM029943.2", "CM029944.2", "CM033063.1", "CM033064.1", "CM029945.2", "CM029946.2"), labels = c("Chr 2L", "Chr 2R", "Chr 3L", "Chr 3R", "Chr XL", "Chr XR"))

# Plot G' analysis without multiple correction p values
gprimeplot <- plotQTLStats(gprime, var = "Gprime", plotThreshold = T) + theme_minimal() 

gprimedist1 <- plotGprimeDist(SNPset = gprime, outlierFilter = "Hampel")  +
  theme_minimal()
gprimedist2 <- plotGprimeDist(SNPset = gprime, outlierFilter = "deltaSNP", filterThreshold = 0.05)  +
  theme_minimal()



# Plot QTL statistics
qtl_offspring <- plotQTLStats(SNPset = qtl_off, var = "deltaSNP", plotIntervals = TRUE) +
  theme_minimal() + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values = c("coral2", "blue"))

##### Trial for better QTL plot- you can choose which one to use

# Plot QTL statistics
qtl_offspring <- plotQTLStats(SNPset = qtl_off, var = "deltaSNP", plotIntervals = TRUE) +
  theme_minimal() + 
  theme(
    axis.title.x = element_text(size = 12),  # Make X-axis title bigger
    axis.title.y = element_text(size = 12),  # Make Y-axis title bigger
    legend.position = "bottom",  # Move legend to the bottom
    legend.title = element_text(face = "bold"),  # Make legend title bold
    legend.text = element_text(face = "bold"),   # Make legend text bold
    strip.text = element_text(face = "bold", size = 10)  # Make facet (chromosome) titles bold and bigger
    
  ) + 
  scale_color_manual(values = c("coral2", "blue"))

# Save the plots

qtl_offspring
ggsave("QTL_map_CI.pdf", plot = qtl_offspring, device = "pdf", width = 20, height = 10)

SNP_density_plot1 <- plotQTLStats(SNPset = qtl_off, var = "nSNPs") + theme_minimal()
ggsave("SNPdensityplot.pdf", plot = SNP_density_plot1, device = "pdf", width = 10, height = 8 )


# Identify significant regions using Delta SNP index and G' methods
sigRegions_qtl <- getSigRegions(qtl_off, method = "QTLseq")
sigRegions_gprime <- getSigRegions(gprime, method = "Gprime")

# Output the number of significant regions for each test
paste0(length(rownames(summary(sigRegions_qtl))), " significant QTL regions")
paste0(length(rownames(summary(sigRegions_gprime))), " significant G' regions")


