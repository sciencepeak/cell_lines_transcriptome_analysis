rm(list = ls())
ptm <- proc.time()
# proc.time() - ptm
options(stringsAsFactors = F)
library("magrittr")

terra_parameter_directory <-
    file.path("..", "metadata") %T>%
    dir.create(., showWarnings = F, recursive = T)

output_table_path <-
    file.path(terra_parameter_directory, "workspace_table.tsv")


# input is a tab delimited table.
raw_file_text <- "
    hisat_index_file    gs://reference_data_02/hisat2/grch38_snp_tran/genome_index.txt
    hisat_prefix	reference_data_02/hisat2/grch38_snp_tran/genome_snp_tran
    ref_gtf gs://reference_data_02/ensembl/Homo_sapiens.GRCh38.99.gtf
    remove_unmapped_inline_single_end_script    gs://reference_data_02/scripts/RemoveUnmappedInlineSingleEnd.pl
    remove_unmapped_inline_paired_end_script    gs://reference_data_02/scripts/RemoveUnmappedInlinePairedEnd.pl
    sorted_bam_result_directory    gs://fc-secure-b93dcb7d-39d9-47ff-9bd9-b281cb9a69de/results/sorted_bam_result_directory
    htseq_count_result_directory  gs://fc-secure-b93dcb7d-39d9-47ff-9bd9-b281cb9a69de/results/htseq_count_result_directory
"

preliminary_workspace_table_dataframe <- read.table(text = raw_file_text)
workspace_table_dataframe <- t(preliminary_workspace_table_dataframe)
# The first column must follow the Terra style:
workspace_table_dataframe[1, 1] <- paste("workspace", workspace_table_dataframe[1, 1], sep = ":")

# -------------------------------------------------------------
# Write hisat_index_file to the hard disk.
# The written file contains address or path to the index files.
# The end of line of the file that will be read by WDL+CromWell+Terra must have 
# Linux style ending, namely "LF". If the Windows style end of line (CRLF) is used,
# The regex error will pop up.
Linux_sytle_output_file <- file(output_table_path, open="wb")

write.table(workspace_table_dataframe,
            file = Linux_sytle_output_file,
            sep = "\t", quote = F,
            col.names = F, row.names = F, eol = "\n")

close(Linux_sytle_output_file)


proc.time() - ptm