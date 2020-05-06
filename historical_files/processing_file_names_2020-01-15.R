rm(list = ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("magrittr")

# -------------------------------------------------------------------------
# --------------------workflow related variables---------------------------

fastq_text <- "M229-DDR-1_S10_R1_001.fastq.gz
M229-DDR-1_S10_R2_001.fastq.gz
M229-DTPPy-2_S1_R1_001.fastq.gz
M229-DTPPy-2_S1_R2_001.fastq.gz
M229-MRDDR2-1_S4_R1_001.fastq.gz
M229-MRDDR2-1_S4_R2_001.fastq.gz
M229-MRDDR4-1_S2_R1_001.fastq.gz
M229-MRDDR4-1_S2_R2_001.fastq.gz
M238-CTRL-1hr_S23_R1_001.fastq.gz
M238-CTRL-1hr_S23_R2_001.fastq.gz
M238-DTPP-2_S5_R1_001.fastq.gz
M238-DTPP-2_S5_R2_001.fastq.gz
M238-DTPPy-2_S6_R1_001.fastq.gz
M238-DTPPy-2_S6_R2_001.fastq.gz
M395-RDDR-1_S9_R1_001.fastq.gz
M395-RDDR-1_S9_R2_001.fastq.gz
M397-Ctrl-4-Day_S13_R1_001.fastq.gz
M397-Ctrl-4-Day_S13_R2_001.fastq.gz
M397-DTPP2-2_S8_R1_001.fastq.gz
M397-DTPP2-2_S8_R2_001.fastq.gz
M397-DTPPy_S24_R1_001.fastq.gz
M397-DTPPy_S24_R2_001.fastq.gz
M397-IFNy-4-Day_S18_R1_001.fastq.gz
M397-IFNy-4-Day_S18_R2_001.fastq.gz
M397-PLX-AZD-11-Day_S16_R1_001.fastq.gz
M397-PLX-AZD-11-Day_S16_R2_001.fastq.gz
M397-PLX-AZD-14-Day_S17_R1_001.fastq.gz
M397-PLX-AZD-14-Day_S17_R2_001.fastq.gz
M397-PLX-AZD-4-Day_S14_R1_001.fastq.gz
M397-PLX-AZD-4-Day_S14_R2_001.fastq.gz
M397-PLX-AZD-7-Day_S15_R1_001.fastq.gz
M397-PLX-AZD-7-Day_S15_R2_001.fastq.gz
M397-PLX-AZD-IFNy-11-Day_S21_R1_001.fastq.gz
M397-PLX-AZD-IFNy-11-Day_S21_R2_001.fastq.gz
M397-PLX-AZD-IFNy-14-Day_S22_R1_001.fastq.gz
M397-PLX-AZD-IFNy-14-Day_S22_R2_001.fastq.gz
M397-PLX-AZD-IFNy-4-Day_S19_R1_001.fastq.gz
M397-PLX-AZD-IFNy-4-Day_S19_R2_001.fastq.gz
M397-PLX-AZD-IFNy-7-Day_S20_R1_001.fastq.gz
M397-PLX-AZD-IFNy-7-Day_S20_R2_001.fastq.gz
SKMEL28-DDR-1_S11_R1_001.fastq.gz
SKMEL28-DDR-1_S11_R2_001.fastq.gz
SKMEL28-DTPPy1-1_S3_R1_001.fastq.gz
SKMEL28-DTPPy1-1_S3_R2_001.fastq.gz
SKMEL28-MRDDR-1_S12_R1_001.fastq.gz
SKMEL28-MRDDR-1_S12_R2_001.fastq.gz
SKMEL28-MRDDR2-1_S7_R1_001.fastq.gz
SKMEL28-MRDDR2-1_S7_R2_001.fastq.gz"


index_file_text <- "genome_snp_tran.1.ht2
genome_snp_tran.2.ht2
genome_snp_tran.3.ht2
genome_snp_tran.4.ht2
genome_snp_tran.5.ht2
genome_snp_tran.6.ht2
genome_snp_tran.7.ht2
genome_snp_tran.8.ht2
"

hisat2_index_folder <- "grch38_snp_tran"
hisat2_index_direct_prefix <- "genome_snp_tran"
genome_index_file_name <- "genome_index.txt"

bucket_address <- "gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c"
bucket_name <- "fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c"

# --------------------workflow related variables---------------------------
# -------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Read in the fastq file and index file names for later aggregation and concatenation.
fastq_file_names <- read.table(text = fastq_text) %>%
    unlist %>%
    unname

index_file_names <- read.table(text = index_file_text) %>%
    unlist %>%
    unname
# ------------------------------------------------------------------------------


# Get the sample names from the fastq files.
sample_names <- sapply(fastq_file_names, function(x) strsplit(x, "_")[[1]][1]) %>%
    unname

# Concatenate the address and fastq file names.
fastq_file_paths <- file.path(bucket_address, "samples", fastq_file_names)

# Group the fastq file paths (google cloud bucket url) based on their sample names
# So we have address for all the r1 fastq files and r2 fastq files for each sample name.
r1_r2_fastq_paths <- tapply(fastq_file_paths, INDEX = sample_names, function(x) paste(x, collapse = ";"))
r1_fastq_paths <- sapply(r1_r2_fastq_paths, function(x) strsplit(x, ";")[[1]][1])
r2_fastq_paths <- sapply(r1_r2_fastq_paths, function(x) strsplit(x, ";")[[1]][2])


# -------------------------------------------------------------
# Write hisat_index_file to the hard disk.
# The written file contains address or path to the index files.
# The end of line of the file that will be read by WDL+CromWell+Terra must have 
# Linux style ending, namely "LF". If the Windows style end of line (CRLF) is used,
# The regex error will pop up.
Linux_sytle_output_file <- file(genome_index_file_name, open="wb") 
output_index_paths <- file.path(bucket_address, hisat2_index_folder, index_file_names)
write.table(output_index_paths, file = Linux_sytle_output_file, quote = F, col.names = F, row.names = F, eol = "\n")
close(Linux_sytle_output_file)


# Assemble all the trivial column columns to fit in with the data table on Terra WorkSpace.
# The column names follow the sytle of Terra and Katie

entity_sample_id <- unique(sample_names)

base_file_name <- unique(sample_names)

fastqTsvFile <- "" %>%
    rep(., times = length(unique(sample_names)))

hisat_index_file <- file.path(bucket_address, hisat2_index_folder, genome_index_file_name) %>%
    rep(., times = length(unique(sample_names)))

hisat_prefix <- file.path(bucket_name, hisat2_index_folder, hisat2_index_direct_prefix) %>%
    rep(., times = length(unique(sample_names)))

input_type <- "fastq" %>%
    rep(., times = length(unique(sample_names)))

# Why is this gtf file needed?
reference_gtf <- "gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/Homo_sapiens.GRCh38.94.gtf" %>%
    rep(., times = length(unique(sample_names)))

# Our library preparation all follows "RF" strandness
strandness <- "RF" %>%
    rep(., times = length(unique(sample_names)))

fastqr1 <- r1_fastq_paths

fastqr2 <-  r2_fastq_paths

all_meta_data_df <- cbind(entity_sample_id,
                          base_file_name,
                          fastqTsvFile,
                          hisat_index_file,
                          hisat_prefix,
                          input_type,
                          reference_gtf,
                          strandness,
                          fastqr1,
                          fastqr2)

# The first column must follow the Terra style: entity:xxx_id
colnames(all_meta_data_df)[colnames(all_meta_data_df) == "entity_sample_id"] <- "entity:sample_id"

# This table needs to be uploaded to Terra/WorkSpace/DATA/Tables on the website.
write.table(all_meta_data_df, file = "hisat_metadata_table.tsv", quote = F, sep = "\t", row.names = F)


proc.time() - ptm


# -----------------------Garbage code:------------------------------------

# # fastqTsvFile_txt_path
# gs_fastqTsvFile_txt_paths <- paste(unique(sample_names), "txt", sep = ".") %>%
#     file.path(bucket_address, "samples", "fastqTsvFile", .)
# 
# # Write fastqTsvFile files to hard disk as Katie did.
# fastqTsvFile_path <- file.path(".", "fastqTsvFile")
# 
# if (dir.exists(fastqTsvFile_path)){
#     unlink(fastqTsvFile_path, recursive = TRUE)
# }
# dir.create(fastqTsvFile_path, showWarnings = FALSE)
# 
# 
# for (i in seq_along(unique(sample_names))) {
#     
#     output_file_name <- unique(sample_names)[i] %>%
#         paste(., "txt", sep = ".")
#     output_path_to_file <- file.path(".", "fastqTsvFile", output_file_name)
#     
#     output_fastqr1_string <- r1_fastq_paths[i]
#     output_fastqr2_string <- r2_fastq_paths[i]
#     cat(output_fastqr1_string, output_fastqr2_string, file = output_path_to_file, sep = "\n")
# }

# index_file_text <- "genotype_genome.1.ht2
# genotype_genome.2.ht2
# genotype_genome.3.ht2
# genotype_genome.4.ht2
# genotype_genome.5.ht2
# genotype_genome.6.ht2
# genotype_genome.7.ht2
# genotype_genome.8.ht2
# genotype_genome.allele
# genotype_genome.clnsig
# genotype_genome.coord
# genotype_genome.fa
# genotype_genome.fa.fai
# genotype_genome.haplotype
# genotype_genome.index.snp
# genotype_genome.link
# genotype_genome.locus
# genotype_genome.partial
# genotype_genome.snp"
# 

# # Group the file names based on their sample names
# r1_r2_fastq_files <- tapply(fastq_file_names, INDEX = sample_names, function(x) paste(x, collapse = ";"))
# r1_fastq_files <- sapply(r1_r2_fastq_files, function(x) strsplit(x, ";")[[1]][1])
# r2_fastq_files <- sapply(r1_r2_fastq_files, function(x) strsplit(x, ";")[[1]][2])
