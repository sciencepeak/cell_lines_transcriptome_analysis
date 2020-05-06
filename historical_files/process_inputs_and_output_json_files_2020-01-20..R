rm(list = ls())
ptm <- proc.time()
# proc.time() - ptm
options(stringsAsFactors = F)
library("magrittr")

# install.packages("jsonlite")
library("jsonlite")

# setwd("D:\\MyData\\Git\\RNAseq_Pipeline\\scripts")


# This will create a file called myWorkflow_inputs.json
# the json file lists all the inputs to all the tasks in your script following the pattern below:
# {"<workflow name>.<task name>.<variable name>": "<variable type>"}

# In the case of Terra, the json file looks like this:
# {"SimplifiedHisatAlignmentWorkflow.strandness":"${this.strandness}",
# "SimplifiedHisatAlignmentWorkflow.fastqr2":"${this.fastqr2}",
# "SimplifiedHisatAlignmentWorkflow.ref_pac":"${workspace.ref_pac}"

# We need to distinguish the entity variables and workspace variables.
# The entity variable has a prefix: this
# The workspace variable has a prefix: workspace


# The user needs to input the following variables and distinguish between entity and workspace variables.
# -----------------------------------------------------------------------
# ---------------------------User write the following names---------------
workflow_name <- "SimplifiedHisatAlignmentWorkflow"

inputs_entity_text <- "
    sample_id
    base_file_name
    strandness
    fastqr1
    fastqr2
    "

inputs_workspace_text <- "
    hisat_index_file
    hisat_prefix
    ref_dict
    ref_fasta
    ref_amb
    ref_ann
    ref_bwt
    ref_fasta_index
    ref_pac
    ref_sa
    dbsnp
    dbsnp_index
    known_indels
    known_indels_index
    thousG
    thousG_index
    "
output_entity_text <- "
    finalBam
    mapped_sam
    md_bam
    sorted_bam
    "


# ---------------------------User write the above names------------------
# -----------------------------------------------------------------------


PasteToList <- function(terra_variables_text, metadata_type) {
    
    metadata_variables <- read.table(text = terra_variables_text) %>%
        unlist %>%
        unname
    
    if (metadata_type == "entity") {
        metadata_prefix <- "this"
    } else if (metadata_type == "workspace") {
        metadata_prefix <- "workspace"
    } else {
        stop("The script can only process entity or workspace metadata")
    }
    
    metadata_values <- paste0("${", metadata_prefix, ".", metadata_variables, "}")
    metadata_keys <- paste(workflow_name, metadata_variables, sep = ".")
    encoded_list <- as.list(metadata_values)
    names(encoded_list) <- metadata_keys
    
    return(encoded_list)
}

inputs_entity_list <- PasteToList(inputs_entity_text, "entity")
inputs_workspace_list <- PasteToList(inputs_workspace_text, "workspace")
outputs_entity_list <- PasteToList(output_entity_text, "entity")

inputs_merged_list <- append(inputs_entity_list, inputs_workspace_list)
outputs_merged_list <- outputs_entity_list

inputs_merged_json <- toJSON(inputs_merged_list, auto_unbox = T)
outputs_merged_json <- toJSON(outputs_merged_list, auto_unbox = T)

cat(inputs_merged_json, file = "../metadata/inputs___.json")
cat(outputs_merged_json, file = "../metadata/outputs___.json")



# In case you want to read a json file
# resource_text <- "../metadata/inputs.json"
# resource_list <- fromJSON(resource_text)









proc.time() - ptm
