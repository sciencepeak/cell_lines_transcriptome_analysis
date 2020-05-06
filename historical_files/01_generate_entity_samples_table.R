rm(list = ls())
ptm <- proc.time()
# proc.time() - ptm
options(stringsAsFactors = F)
library("magrittr")

# --------------------workflow related variables---------------------------

bucket_address <- "gs://fc-secure-9b93eac8-660e-4bf4-a103-7d96d7054601"

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


# ---------------------------------------------------------------------------
# Read in the fastq file and index file names for later aggregation and concatenation.
fastq_file_names <- read.table(text = fastq_text) %>%
    unlist %>%
    unname

# Get the sample names from the fastq files.
sample_names <- sapply(fastq_file_names, function(x) strsplit(x, "_")[[1]][1]) %>%
    unname

# Concatenate the address and fastq file names.
fastq_file_paths <- file.path(bucket_address, "samples", fastq_file_names)


# ---------------------------------------------------------------------------------
# Group the fastq file paths (google cloud bucket url) based on their sample names
# So we have address for all the r1 fastq files and r2 fastq files for each sample name.
r1_r2_fastq_paths <- tapply(fastq_file_paths, INDEX = sample_names, function(x) paste(x, collapse = ";"))
r1_fastq_paths <- sapply(r1_r2_fastq_paths, function(x) strsplit(x, ";")[[1]][1])
r2_fastq_paths <- sapply(r1_r2_fastq_paths, function(x) strsplit(x, ";")[[1]][2])

names(r1_r2_fastq_paths) == unique(sample_names)
names(r1_fastq_paths) == unique(sample_names)
names(r2_fastq_paths) == unique(sample_names)
sample_names == sort(sample_names)
unique(sample_names) == sort(unique(sample_names))

# Assemble all the trivial column columns to fit in with the data table on Terra WorkSpace.
# The column names follow the sytle of Terra and Katie

entity_sample_id <- unique(sample_names)

base_file_name <- unique(sample_names)

# Most pair-end library preparation follows "RF" strandness
# Our data is unstranded.
# strandness <- rep("unstranded", times = length(sample_names))

fastqr1 <- r1_fastq_paths

fastqr2 <-  r2_fastq_paths

all_meta_data_df <- cbind(entity_sample_id,
                          base_file_name,
                          fastqr1,
                          fastqr2)

# The first column must follow the Terra style: entity:xxx_id
colnames(all_meta_data_df)[colnames(all_meta_data_df) == "entity_sample_id"] <- "entity:sample_id"

# This table needs to be uploaded to Terra/WorkSpace/DATA/Tables on the website.
write.table(all_meta_data_df, file = "../metadata/entity_samples_table.tsv", quote = F, sep = "\t", row.names = F)


proc.time() - ptm
