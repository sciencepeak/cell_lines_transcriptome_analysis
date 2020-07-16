rm(list = ls())
ptm <- proc.time()
# proc.time() - ptm
options(stringsAsFactors = F)
library("magrittr")
library("tidyr")
library("stringr")
library("readxl")
library("tibble")


output_file_path <- "../metadata/entity_samples_table.tsv"

# List all files in the google bucket.
batch_01_files <- system2("gsutil", args =  c("ls", "gs://reference_data_02/samples/cell_lines_batch_01"), stdout = TRUE)
batch_02_files <- system2("gsutil", args =  c("ls", "gs://reference_data_02/samples/cell_lines_batch_02"), stdout = TRUE)

# Subset the fastq.gz files.
batch_01_fastq_files <-  subset(batch_01_files, grepl("fastq.gz", batch_01_files))
batch_02_fastq_files <-  subset(batch_02_files, grepl("fastq.gz", batch_02_files))

# Download the annotation grouping file from google cloud bucket.
grouping_table_path <- subset(batch_01_files, grepl("xlsx", batch_01_files))
system2("gsutil", args = c("cp", grouping_table_path, "../metadata"))

# Read the grouping information into R workspace
grouping_table_file <- basename(grouping_table_path)
grouping_dataframe <- file.path("..", "metadata", grouping_table_file) %>%
    read_excel %>%
    as.data.frame %>%
    inset(., "Batch", value = paste("batch", .$Batch, sep = "_"))
    

# Ensure there is no duplicates of sample names in the two batches.
batch_01_fastq_files %in% batch_02_fastq_files
batch_02_fastq_files %in% batch_01_fastq_files

# Combine the two batches of fastq file paths on google bucket.
fastq_file_paths <- c(batch_01_fastq_files, batch_02_fastq_files)

# Extract sample names from the fastq file paths.
# The sample names now have duplicates.
sample_names <- fastq_file_paths %>%
    basename %>%
    str_split(., pattern = "_", simplify = T) %>%
    .[, 1]

# Ensure annotation grouping table's sample names comply with the sample names extracted from fastq files.
grouping_dataframe[, 1] %in% unique(sample_names)

# Group the fastq file paths (google cloud bucket url) based on their sample names
# So we have address for all the r1 fastq files and r2 fastq files for each sample name.
r1_r2_fastq_paths <- tapply(fastq_file_paths,
                            INDEX = sample_names,
                            function(x) paste(x, collapse = ";"))


separate_R1R2_fastq_path <- function(path_mixture_string) {
    
    mixed_dataframe <- path_mixture_string %>%
        str_split(., pattern = ";", simplify = T)
    
    upper_limit <- ncol(mixed_dataframe)
    
    odd_index <- seq(from = 1, to = upper_limit, by = 2)
    even_index <- seq(from = 2, to = upper_limit, by = 2)
    
    odd_paths <- mixed_dataframe[1, odd_index]
    even_paths <- mixed_dataframe[1, even_index]
    
    combined_odd_path <- paste(odd_paths, collapse = ";")
    combined_even_path <- paste(even_paths, collapse = ";")
    
    R1R2_path_vector <- c(combined_odd_path, combined_even_path)
    names(R1R2_path_vector) <- c("fastq_r1", "fastq_r2")
    
    return(R1R2_path_vector)
}

# Concatenate the fastq file paths into R1 and R2.
sample_path_dataframe <- sapply(r1_r2_fastq_paths, separate_R1R2_fastq_path) %>%
    t %>%
    as.data.frame %>%
    rownames_to_column(., "base_file_name")

# Merge the sample and paths with annotation of grouping information and batch information.
all_metadata_dataframe <- merge(sample_path_dataframe, grouping_dataframe, by.x = "base_file_name", by.y = "File", all.x = T) %>%
    unite(., col = "entity:sample_id", c("Batch", "Group", "base_file_name"), sep = "_", remove = F) %$%
    .[order(Batch, base_file_name, Group), ] %>%
    .[, 1:4]

# Write out entity table to the hard disk.
# This table needs to be uploaded to Terra/WorkSpace/DATA/Tables on the website.
write.table(all_metadata_dataframe,
            file = output_file_path,
            quote = F, sep = "\t", row.names = F)


proc.time() - ptm
