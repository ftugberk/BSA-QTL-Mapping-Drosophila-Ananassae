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

### Annotation Preparation

# Import gene annotation data and ortholog information
refseq_flybase <- read_tsv("REFSEQ_FLYBASE_Dana.txt") # contains REFSEQ IDs and FLYBASE ID
dmel_dana_ortho <- read_excel("dmel_dana_orthologs.xlsx") # contains D.melanogaster orthologs of D.ananassae

# Select relevant columns and join with ortholog data
orthologs <- AnnotationDbi::select(org.Dm.eg.db, keys = dmel_dana_ortho$ortholog, keytype = "FLYBASE", columns = c("FLYBASE", "ENTREZID")) %>%
  right_join(dmel_dana_ortho, by = c("FLYBASE" = "ortholog"))


annotation <- import("genomic.gtf") # NCBI Drosophila ananassae REFSEQ Annotation Release 102

# Get significant QTL regions data with getQTLTable function from QTLseqr
sigQTL <- getQTLTable(qtl_off, method = "QTLseq", export = T, fileName = "sigQTL.csv")

# Mapping of chromosome names
chromosome_mapping <- c(
  "Chr 2L" = "NC_057927.1",
  "Chr 2R" = "NC_057928.1",
  "Chr 3L" = "NC_057929.1",
  "Chr 3R" = "NC_057930.1",
  "Chr XL" = "NC_057931.1",
  "Chr XR" = "NC_057932.1"
)


# Change chromosome names to their originals
sigQTL <- sigQTL %>%
  mutate(CHROM = recode(CHROM, !!!chromosome_mapping))

# Divide the QTL regions into positive and negative based on Delta SNP index
pos_QTL <- sigQTL[which(sigQTL$avgDeltaSNP > 0),1:4]
neg_QTL <- sigQTL[which(sigQTL$avgDeltaSNP < 0),1:4]


# Function to process QTL data and find overlapping genes
process_qtl <- function(sigQTL_subset, annotation, refseq_flybase, orthologs, output_file) {
  
  # Convert QTL subset to GRanges object
  granges_sig_regions <- makeGRangesFromDataFrame(sigQTL_subset[,-2], keep.extra.columns = TRUE)
  
  # Find overlaps with the annotation data
  overlaps <- findOverlaps(granges_sig_regions, annotation)
  overlapping_genes <- annotation[subjectHits(overlaps)]
  
  # Extract unique transcript IDs
  unique_transcript_ids <-
    as.data.frame(unique(overlapping_genes@elementMetadata@listData[["transcript_id"]]))
  colnames(unique_transcript_ids) <- "transcript_id"
  
  # Split the transcript_id into REFSEQ components and filter for relevant REFSEQ types
  split_transcript_ids <- unique_transcript_ids %>%
    separate(transcript_id, c("REFSEQ_TYPE", "REFSEQ_ID"), "_") %>%
    filter(REFSEQ_TYPE == "XM")
  
  # Reunite the REFSEQ_TYPE and REFSEQ_ID into a single column
  filtered_transcripts <- split_transcript_ids %>%
    unite("REFSEQ", REFSEQ_TYPE, REFSEQ_ID, sep = "_")
  
  # Remove version numbers from REFSEQ in filtered_transcripts
  filtered_transcripts$REFSEQ <- sub("\\.\\d+$", "", filtered_transcripts$REFSEQ)
  
  # Merge with refseq_flybase data based on REFSEQ to add additional identifiers (Dana_ID)
  annotated_transcripts <- merge(x = filtered_transcripts, y = refseq_flybase, by = "REFSEQ")
  colnames(annotated_transcripts) <- c("REFSEQ", "Dana_ID")
  
  # Merge with ortholog data to map Dana_ID to Drosophila melanogaster orthologs
  annotated_orthologs <- merge(x = annotated_transcripts, y = orthologs, by = "Dana_ID")
  
  # Write the final annotated orthologs to a TSV file
  write_tsv(annotated_orthologs, output_file)
}


# Process positive, negative, and all QTLs
ogenes_dmelpositive <- process_qtl(pos_QTL, annotation, refseq_flybase, orthologs, "positive_annotated_orthologs.tsv")
ogenes_dmelnegative <- process_qtl(neg_QTL, annotation, refseq_flybase, orthologs, "negative_annotated_orthologs.tsv")

# Load necessary libraries for this part

library(ggplot2)
library(clusterProfiler)
library(org.Dm.eg.db)  

# Define a function to perform GO enrichment and plot results
perform_GO_analysis_tubi <- function(result_data, org_db, region_label) {
  # Perform GO enrichment for BP, MF, and CC
  gseGO_BP <- perform_enrichGO_tubi(result_data$ENTREZID, "BP", org_db)
  gseGO_MF <- perform_enrichGO_tubi(result_data$ENTREZID, "MF", org_db)
  gseGO_CC <- perform_enrichGO_tubi(result_data$ENTREZID, "CC", org_db)
  
  # Process the results
  BP <- process_enrich_result_tubi(gseGO_BP, "Biological Process", region_label)
  MF <- process_enrich_result_tubi(gseGO_MF, "Molecular Function", region_label)
  CC <- process_enrich_result_tubi(gseGO_CC, "Cellular Component", region_label)
  
  # Combine data
  combined_data <- rbind(BP, MF, CC)
  
  return(combined_data)
}

# Define a function to perform GO enrichment
perform_enrichGO_tubi <- function(entrez_ids, ont, org_db) {
  enrichGO(gene = entrez_ids,
           OrgDb = org_db,
           keyType = "ENTREZID",
           ont = ont,
           readable = T,
           pAdjustMethod = "fdr",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)
}

# Define a function to process enrichment results
process_enrich_result_tubi <- function(enrich_result, ontology, region_label, p_adjust_cutoff = 0.05) {
  result_df <- as.data.frame(enrich_result@result)
  result_df <- result_df[result_df$p.adjust <= p_adjust_cutoff, ]
  result_df$Ontology <- ontology
  result_df$Region <- region_label
  return(result_df)
}

# Positive and negative regions GO analysis
positive_data_tubi <- perform_GO_analysis_tubi(ogenes_dmelpositive, org.Dm.eg.db, "Positive")
negative_data_tubi <- perform_GO_analysis_tubi(ogenes_dmelnegative, org.Dm.eg.db, "Negative")

# Combine positive and negative data for combined plot
both_data_tubi <- rbind(positive_data_tubi, negative_data_tubi)


### Plots

# Plot for Negative Region
negative_plot_tubi <- ggplot(negative_data_tubi, aes(x = reorder(Description, -Count), y = Count, fill = Region)) +
  geom_bar(stat = "identity", fill = "tomato", position = "dodge") +
  facet_grid(. ~ Ontology, scales = "free", space = "free") +
  scale_y_continuous(breaks = seq(0, 80, by = 5)) + 
  labs(title = "Negative QTLs GO Enrichment Analysis", x = "GO Terms", y = "Number of Genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(color = "black", fill = "white", size = 1.5, linetype = "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

# Plot for Positive Region
positive_plot_tubi <- ggplot(positive_data_tubi, aes(x = reorder(Description, -Count), y = Count, fill = Region)) +
  geom_bar(stat = "identity", fill = "steelblue", position = "dodge") +
  facet_grid(. ~ Ontology, scales = "free", space = "free") +
  scale_y_continuous(breaks = seq(0, 80, by = 5)) + 
  labs(title = "Positive QTLs GO Enrichment Analysis", x = "GO Terms", y = "Number of Genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(color = "black", fill = "white", size = 1.5, linetype = "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20),
    axis.title = element_text(face = "bold",size = 12),
    legend.position = "none"
  )

# Plot for both regions
combined_plot_tubi <- ggplot(both_data_tubi, aes(x = reorder(Description, -Count), y = Count, fill = Region)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Ontology, scales = "free", space = "free") +
  scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +
  scale_y_continuous(breaks = seq(0, 80, by = 5)) + 
  labs(title = "Gene Ontology Enrichment Analysis", x = "Gene Ontology Terms", y = "Number of Genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),  
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(color = "black", fill = "white", size = 1.5, linetype = "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 50),  
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  )



# Save the plots
ggsave("negative_GO_plot_tubi.pdf", plot = negative_plot_tubi, device = "pdf", width = 15, height = 10)
ggsave("positive_GO_plot_tubi.pdf", plot = positive_plot_tubi, device = "pdf", width = 15, height = 10)
ggsave("combined_GO_plot_tubi.pdf", plot = combined_plot_tubi, device = "pdf", width = 20, height = 10)

# Display the plots
print(negative_plot_tubi)
print(positive_plot_tubi)
print(combined_plot_tubi)

### Supplementary Dot Plots

#Convert 'GeneRatio' to numeric by evaluating the fraction, then add 'GeneRatioPercent' 

positive_data_tubi <- positive_data_tubi %>%
  mutate(GeneRatio = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
         GeneRatioPercent = GeneRatio * 100)  # Convert to percentage

negative_data_tubi <- negative_data_tubi %>%
  mutate(GeneRatio = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
         GeneRatioPercent = GeneRatio * 100)  # Convert to percentage

# Convert Description to factor with reversed levels for better plotting
positive_data_tubi$Description <- factor(positive_data_tubi$Description, levels = rev(unique(positive_data_tubi$Description)))
negative_data_tubi$Description <- factor(negative_data_tubi$Description, levels = rev(unique(negative_data_tubi$Description)))

#Create Dot Plot for Positive Region

dot_plot_positive <- ggplot(positive_data_tubi, aes(x = GeneRatioPercent, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) + 
  scale_size_continuous(name = "Number of Genes") +  # Size based on gene count
  scale_color_gradient(name = "Adjusted p-value", low = "red", high = "blue") +  # Color gradient for p-value
  theme_minimal() + 
  labs(title = "Gene Ontology Enrichment for Positive QTL", x = "Gene Ratio (%)", y = "Gene Ontology Terms") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  )

# Display the positive plot
print(dot_plot_positive)

# Save the positive plot
ggsave("positive_go_enrichment_dot_plot.png", plot = dot_plot_positive, width = 12, height = 8)

# Step 3: Create Dot Plot for Negative Region
dot_plot_negative <- ggplot(negative_data_tubi, aes(x = GeneRatioPercent, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) + 
  scale_size_continuous(name = "Number of Genes") +  # Size based on gene count
  scale_color_gradient(name = "Adjusted p-value", low = "red", high = "blue") +  # Color gradient for p-value
  theme_minimal() + 
  labs(title = "Gene Ontology Enrichment for Negative QTL", x = "Gene Ratio (%)", y = "Gene Ontology Terms") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  )

# Display the negative plot
print(dot_plot_negative)

# Save the negative plot
ggsave("negative_go_enrichment_dot_plot.png", plot = dot_plot_negative, width = 12, height = 8)

