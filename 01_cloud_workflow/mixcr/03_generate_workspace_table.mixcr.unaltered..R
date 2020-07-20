rm(list = ls())
ptm <- proc.time()
# proc.time() - ptm
options(stringsAsFactors = F)
library("magrittr")

terra_parameter_directory <-
    file.path("..", "terra_parameters.mixcr") %T>%
    dir.create(., showWarnings = F, recursive = T)

output_table_path <-
    file.path(terra_parameter_directory, "workspace_table.mixcr.unaltered.tsv")

raw_file_text <- "
    mixcr_clones_result_directory   gs://fc-secure-b93dcb7d-39d9-47ff-9bd9-b281cb9a69de/results/mixcr_clones_result_directory
    mixcr_suffix    mixcr_clones.txt
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