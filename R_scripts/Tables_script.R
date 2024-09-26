### This script is to creat necessary supplementary Tables and Excel Sheets

### Excel sheet on GO terms and Gene IDs

# Load library

library(writexl)
library(readr)
library(kableExtra)
library(knitr)



# Filter the data frame to keep only 'Description' and 'geneID' columns
GeneID_combined_data <- both_data_tubi[, c("Description", "geneID")]

# Save the filtered data frame to an Excel file
write_xlsx(GeneID_combined_data, "GeneID_combine_data.xlsx")

### Table: Significant QTL tables



# Step 2: Create a lookup table for chromosome names
chromosome_lookup <- c(
  "NC_057927.1" = "Chr 2L",
  "NC_057928.1" = "Chr 2R",
  "NC_057929.1" = "Chr 3L",
  "NC_057930.1" = "Chr 3R",
  "NC_057931.1" = "Chr XL",
  "NC_057932.1" = "Chr XR"
)

# Step 3: Replace the CHROM column values using the lookup table
sigQTL$CHROM <- chromosome_lookup[sigQTL$CHROM]

# Step 4: Select relevant columns
qtl_table <- sigQTL[, c("CHROM", "start", "end", "length", "nSNPs", "avgDeltaSNP")]

sig_qtl_table <- qtl_table %>%
  kable("html", caption = "Significant QTL Regions") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                full_width = FALSE, 
                position = "center")

### Table: Variant Counts on vcf file through stages of filtration

VCF_File_Statistics_Combined <- read_csv("C:/Users/ftugb/OneDrive/Desktop/Master Thesis/R_analysis/VCF_File_Statistics_Combined.csv")

# Create a formatted table
table <- kbl(VCF_File_Statistics_Combined, booktabs = TRUE, caption = "VCF File Statistics") %>%
  kable_styling(full_width = FALSE, position = "center", font_size = 12) %>%
  column_spec(1, bold = TRUE) # Make the first column bold

# Save the table as an image or HTML file
save_kable(table, file = "VCF_File_Statistics_Table.html")

### Table: Quality metrics of the raw sequincing data files 


# Create a data frame with the table's data, removing "Raw_data" column
QC_rawdata <- data.frame(
  Sample = c("slow O", "slow O", "slow O", "fast O", "fast O"),
  Library_Flowcell_Lane = c("EKDN230016689 1A HF3MMDSX7 L1", 
                            "EKDN230016689 1A HF37MDSX7 L3", 
                            "EKDN230016689 1A HF35WDSX7 L3", 
                            "EKDN230016690 1A HF3WTDSX7 L2", 
                            "EKDN230016690 1A HF7TMDSX7 L3"),
  Raw_reads = c(14739612, 5796120, 100075174, 89512044, 37727344),
  Effective = c(98.78, 98.83, 98.84, 98.28, 98.43),
  Error = c(0.03, 0.03, 0.03, 0.03, 0.03),
  Q20 = c(95.87, 96.98, 97.46, 96.43, 97.12),
  Q30 = c(89.85, 92.14, 92.98, 91.08, 92.37),
  GC = c(45.14, 45.36, 45.00, 43.67, 43.72)
)

# Create the table using kable
QC_raw_data_table <- kable(QC_rawdata, "html", col.names = c("Sample", "Library Flowcell Lane", "Raw reads", 
                                                             "Effective (%)", "Error (%)", "Q20 (%)", "Q30 (%)", "GC (%)")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  column_spec(1, bold = TRUE) %>%
  add_header_above(c(" " = 1, "Details" = 2, "Quality Metrics" = 5))


