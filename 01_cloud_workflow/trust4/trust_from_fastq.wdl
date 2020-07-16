version 1.0
workflow MyBestWorkflow {
    input {
        # the sample_id in an entity table stores some raw information of the samples
        # The base_file_name is the actual, concise, and actual sample_name throughout the code.
        String base_file_name

        # Required for (;-concatednated) fastq files
        File fastq_r1
        File? fastq_r2

        Boolean? single_endness

        String trust_suffix
        String trust_clones_result_directory
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

    call TrustFromFastq {
        input:
            sample_name = base_file_name,
            fastq_files = ConcatenateGzippedFastq.concatenated_gzipped_fastq_files,
            single_endness = single_endness
    }

    call CompressFileToBzip2 {
        input:
            sample_name = base_file_name,
            input_file = TrustFromFastq.trust_report_txt,
            appended_annotation = trust_suffix
    }

    call CollectCloudFiles {
        input:
            source_file = CompressFileToBzip2.compressed_bzip2_file,
            destination_directory = trust_clones_result_directory
    }

    # Output files of the workflows.
    output {
        File trust_outfile = CompressFileToBzip2.compressed_bzip2_file
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

task TrustFromFastq {
    input {
        String sample_name
        Array[File]+ fastq_files
        Boolean? single_endness
    }

    File fastq_r1_file = fastq_files[0]
    File? fastq_r2_file = fastq_files[1]

    Boolean single_end_argument = select_first([single_endness, false])

    command {
        if [[ "${single_end_argument}" == true ]]
            then
                echo "The single end input is detected"
                /home/TRUST4/run-trust4 -t 1 -u ${fastq_r1_file} -f /home/TRUST4/hg38_bcrtcr.fa --ref /home/TRUST4/human_IMGT+C.fa -o ${sample_name}
            else
                echo "The paired end input is detected"
                /home/TRUST4/run-trust4 -t 1 -1 ${fastq_r1_file} -2 ${fastq_r2_file} -f /home/TRUST4/hg38_bcrtcr.fa --ref /home/TRUST4/human_IMGT+C.fa -o ${sample_name}
        fi
    }

    output {
        File trust_report_txt = glob("*_report.tsv")[0]
    }

    runtime {
        memory: "32G"
        cpu: 1
        disks: "local-disk 375 SSD"
        docker: "jemimalwh/trust4:0.2.0"
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
