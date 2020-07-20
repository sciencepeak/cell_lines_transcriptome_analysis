version 1.0
workflow MyBestWorkflow {
    input {
        # IMPORTANT:
        # This script only handles fastq. fastqList is not handled.

        # FILE INFO
        # the sample_id stores some raw information of the samples
        # The base_file_name is the actual, concise, and actual sample_name throughout the code.
        String base_file_name

        # Required for (concatednated) fastq files
        File fastq_r1
        File? fastq_r2

        String mixcr_suffix
        # The output files will be copied here on Google Drive Bucket or local server
        String mixcr_clones_result_directory
    }


    call SplitFastqNames {
        input:
            sample_name = base_file_name,
            fastq_r1_string = fastq_r1,
            fastq_r2_string = fastq_r2
    }

    call ConcatenateGzippedFastq {
        input:
            sample_name = base_file_name,
            fastq_list_files = SplitFastqNames.fastq_list_files

    }
    call MixcrFromRawFastq {
        input:
            sample_name = base_file_name,
            input_fastq_files = ConcatenateGzippedFastq.concatenated_gzipped_fastq_files
    }

    call CompressFileToBzip2 {
        input:
            sample_name = base_file_name,
            input_file = MixcrFromRawFastq.mixcr_clones_txt,
            appended_annotation = mixcr_suffix
    }
    call CollectCloudFiles {
        input:
            source_file = CompressFileToBzip2.compressed_bzip2_file,
            destination_directory = mixcr_clones_result_directory
    }

    # Output files of the workflows.
    output {
        File mixcr_clones_txt = MixcrFromRawFastq.mixcr_clones_txt
    }
}

task SplitFastqNames {
    input {
        String sample_name
        String fastq_r1_string
        String? fastq_r2_string
    }

    command {
        echo "${fastq_r1_string}" | tr ";" "\n" > ${sample_name}.fastq_r1_list.txt
        echo "${fastq_r2_string}" | tr ";" "\n" > ${sample_name}.fastq_r2_list.txt
    }

    output {
        Array[File] fastq_list_files = glob("*_list.txt")
    }

    runtime {
        disks: "local-disk 375 SSD"
        docker: "ubuntu:latest"
    }
}

task ConcatenateGzippedFastq {
    # Task Input Declaration
    input {
        String sample_name
        Array[File] fastq_list_files
    }

    # Non-Input Declarations
    Array[File] fastq_r1_list = read_lines(fastq_list_files[0])
    Array[File]? fastq_r2_list = read_lines(fastq_list_files[1])


    Int file_length = length(fastq_r1_list)

    command {
        if [[ "${file_length}" == 1 ]]
            then
                echo "The fastq input is one fastq file for each pair"
                cp ${sep=" " fastq_r1_list} ${sample_name}.R1.fastq.gz
                cp ${sep=" " fastq_r2_list} ${sample_name}.R2.fastq.gz
            else
                echo "The fastq input is more than one fastq file for each pair"
                zcat ${sep=" " fastq_r1_list} | gzip -9 -cvf > ${sample_name}.R1.fastq.gz
                zcat ${sep=" " fastq_r2_list} | gzip -9 -cvf > ${sample_name}.R2.fastq.gz
        fi
    }

    output {
        Array[File] concatenated_gzipped_fastq_files = glob("*R*.fastq.gz")
    }

    runtime {
        disks: "local-disk 375 SSD"
        docker: "ubuntu:latest"
    }
}

task MixcrFromRawFastq {
    input {
        String sample_name
        Array[File] input_fastq_files
    }

    File fastq_file_r1 = input_fastq_files[0]
    File fastq_file_r2 = input_fastq_files[1]

    # Typical MiXCR workflow consists of three main processing steps:
    command {
        # step 1: align: align sequencing reads to reference V, D, J and C genes of T- or B- cell receptors
        /mixcr/mixcr align --species HomoSapiens --report ${sample_name}.mixcr_align_report.txt ${fastq_file_r1} ${fastq_file_r2} ${sample_name}.mixcr_aligned.vdjca

        # step 2: assemble: assemble clonotypes using alignments obtained on previous step (in order to extract specific gene regions e.g. CDR3)
        /mixcr/mixcr assemble --report ${sample_name}.mixcr_assemble_report.txt -a ${sample_name}.mixcr_aligned.vdjca ${sample_name}.mixcr_assembled.clna

        # step 3: export: export alignment (exportAlignments) or clones (exportClones) to human-readable text file
        # export alignments from .clna file
        /mixcr/mixcr exportAlignments ${sample_name}.mixcr_assembled.clna ${sample_name}.mixcr_alignments.txt
        # export clones from .clna file
        /mixcr/mixcr exportClones --preset full ${sample_name}.mixcr_assembled.clna ${sample_name}.mixcr_clones.txt
    }

    output {
        File mixcr_clones_txt = "${sample_name}.mixcr_clones.txt"
    }

    # docker docker pull mxhe/mixcr:3.0.12
    runtime {
        disks: "local-disk 750 SSD"
        docker: "mxhe/mixcr:3.0.12"
    }
}

task CompressFileToBzip2 {
    input {
        String sample_name
        File input_file
        String appended_annotation
    }

    command {
        bzip2 -9 -cvf ${input_file} > ${sample_name}.${appended_annotation}.bz2
    }

    output {
        File compressed_bzip2_file = "${sample_name}.${appended_annotation}.bz2"
    }

    runtime {
        memory: "8G"
        cpu: 1
        disks: "local-disk 375 SSD"
        docker: "ubuntu:latest"
    }
}

task CollectCloudFiles {
    input {
        # Source file type is File
        File source_file
        # The destination directory type is String.
        String destination_directory
    }

    command {
        # This is the feasible way of copying files in wdl on Terra
        gsutil cp ${source_file} ${destination_directory}
    }

    runtime {
        disks: "local-disk 375 SSD"
        docker: "google/cloud-sdk:latest"
    }
}
