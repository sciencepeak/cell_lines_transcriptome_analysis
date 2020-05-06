#' Aggregate expression values into a single matrix
#'
#' Given a list of files, construct a data frame of the expression values
#' of genes/transcripts (rows) across samples (columns).
#' @name summarize_expression
#' @param files Character vector specifying file paths containing gene
#' expression output. Default formatting of files is two column format, with
#' unique gene identifiers in column 1 and expression values in column 2.
#' @param dir Character string specifying directory containing all files.
#' @param pattern Character string specifying pattern in file paths (within
#' directory).
#' @param samples List of character strings to use as sample identifiers. If
#' these aren't provided, the columns are labeled Sample_1, Sample_2, etc. in
#' order of the files provided.
#' @param format Character string specifying the program used to generate gene
#' expression output. "HTSeq" and "Stringtie" are currently supported.
#' @param head Boolean value specifying whether gene expression output files
#' have a header or not. If the format is supported, leave blank. Default is
#' FALSE.
#' @param colnames Character vector specifying column names of gene expression
#' output. If the format is supported, leave blank. Default is 'Gene' and
#' 'Expression.'
#' @param idcol Character or integer vector specifying which column contains
#' unique gene/target identifiers. Default is column 1.
#' @param expressioncol Character string or integer specifying the column
#' containing the expression information. Default is column 1 for unspecified
#' or HTSeq format. Default is 'FPKM' for Stringtie format.
#' @param out Character string specifying the object to output, one of "data" or
#' "file", defaults to "data" (see returns).
#' @param outdir Character string specifying output directory, defaults to
#' current working directory.
#' @param outfile Character string specifying output file name prefix. Default
#' is EXPRESSION_MATRIX.txt.
#' @return One of the following, a dataframe containing a gene expression
#' matrix or a file containing the gene expression matrix.
#' @import tidyverse
#' @importFrom utils read.table
#' @export
#'
#'

summarize_expression <- function(files = NULL, dir = getwd(), pattern = NULL,
                                 samples = NULL, format = NULL, head = FALSE,
                                 colnames = c('Gene','Expression'),
                                 idcol = c('Gene'), expressioncol = 'Expression',
                                 out = 'data', outdir = getwd(),
                                 outfile = 'EXPRESSION_MATRIX') {
  
  # If files are not provided, look for files in given directory
  if(is.null(files)){
    files <- list.files(path = dir, pattern = pattern, full.names = T)
    warning('No files were provided. There were ', length(files),
            ' files found in ', dir, ': ', paste(files, collapse = ", "))
  }

  # Check to make sure that input files exist; returns error if any do not
  check_file_exists <- sapply(files, function(xx) {
    file.exists(as.character(xx))
  })
  if(!all(check_file_exists)) {
    dne_files <- files[which(check_file_exists==F)]
    stop("The following files were not found: ",
         paste(dne_files, collapse = ", "))
  }

  # Check to see if sample names are provided; use file names if not
  if(is.null(samples)) {
    names(files) <- paste0('Sample_',1:length(files))
  } else {
    names(files) <- samples
  }

  # Generate formatting details
  if(is.null(format)) {
    df_colnames <- colnames
  } else if(format == 'HTSeq') {
    # If format is HTSeq, columns are 'target_id' (idcol) and 'counts'
    # (expressioncol), with no header
    df_colnames <- c('target_id', 'counts')
    head <- FALSE
    idcol <- c('target_id')
    expressioncol <- 'counts'
  } else if(format == 'Stringtie') {
    # If format is Stringtie, columns are (see next line), with a header.
    # expressioncol can be TPM or FPKM; default
    df_colnames <- c('Gene_ID','Gene_Name','Reference','Strand','Start','End',
                     'Coverage','FPKM','TPM')
    head <- TRUE
    idcol <- c('Gene_ID','Gene_Name','Reference','Strand','Start','End')
    if(expressioncol == 'Expression'){
      expressioncol <- 'FPKM'
    }
  } else {
    df_colnames <- colnames
  }

  samplecol <- "SAMPLE"
  # Read in expression files
  read_files <- lapply(names(files), function(xx) {
    d <- read.table(files[[xx]], sep = "\t", head = head,
                    stringsAsFactors = FALSE,
                    col.names = df_colnames)
    d <- dplyr::mutate(d, SAMPLE = as.character(xx))
    dplyr::select(d, idcol, expressioncol, samplecol)
  })

  long <- Reduce(rbind, read_files)
  output <- tidyr::spread(long, samplecol, expressioncol)

  return(output)
}

