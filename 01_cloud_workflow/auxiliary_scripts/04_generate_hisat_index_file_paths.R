rm(list = ls())
ptm <- proc.time()
# proc.time() - ptm
options(stringsAsFactors = F)
library("magrittr")


index_file_text <- "
	genome_snp_tran.1.ht2
	genome_snp_tran.2.ht2
	genome_snp_tran.3.ht2
	genome_snp_tran.4.ht2
	genome_snp_tran.5.ht2
	genome_snp_tran.6.ht2
	genome_snp_tran.7.ht2
	genome_snp_tran.8.ht2
"
index_file_names <- read.table(text = index_file_text) %>%
    unlist %>%
    unname

output_index_paths <- file.path("gs://reference_data_02",
                                "hisat2",
                                "grch38_snp_tran",
                                index_file_names)
output_index_paths

# -------------------------------------------------------------
# Write hisat_index_file to the hard disk.
# The written file contains address or path to the index files.
# The end of line of the file that will be read by WDL+CromWell+Terra must have 
# Linux style ending, namely "LF". If the Windows style end of line (CRLF) is used,
# The regex error will pop up.
Linux_sytle_output_file <- file("../metadata/genome_index.txt", open="wb")

write.table(output_index_paths, file = Linux_sytle_output_file, quote = F, col.names = F, row.names = F, eol = "\n")
close(Linux_sytle_output_file)

print("Don't forget uploading the latest genome_index file to the Google Cloud Bucket!")

proc.time() - ptm