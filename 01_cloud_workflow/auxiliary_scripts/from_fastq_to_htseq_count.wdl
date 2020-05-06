workflow MyBestWorkflow {
    # IMPORTANT:
    # This script only handles fastq. fastqList is not handled.

    # FILE INFO
    # The base_file_name is actually the sample_name
    String base_file_name

    # Required for fastq
    File? fastq_r1
    File? fastq_r2

    # Handle single_ended or pair_ended?
    Boolean? single_ended
    Boolean single_end_attribute = select_first([single_ended, false])

    # Process input args: usually the library prepation for strandness is "RF".
    String? strandness


    # REFERENCE FILES for HiSAT2
    File hisat_index_file
    Array[String] hisat_index = read_lines(hisat_index_file)
    String hisat_prefix

    # GTF file is needed by HTSeq-count
    File ref_gtf

    # remove unmapped reads from sam files using custom perl scripts.
    File remove_unmapped_inline_single_end_script
    File remove_unmapped_inline_paired_end_script

    # These two locations on Google Drive Bucket are used to collect results,
    # namely, bam files and htseq-count text files.
    String sorted_bam_result_directory
    String htseq_count_result_directory

    call StringToFile {
        input:
            sample_name = base_file_name,
            fastq_r1_string = fastq_r1,
            fastq_r2_string = fastq_r2
    }

    call FastqToSam {
        input:
            fastq_list_files = StringToFile.fastq_list_files,
            single = single_end_attribute,
            sample_name = base_file_name,
            hisat_index = hisat_index,
            hisat_prefix = hisat_prefix,
            strandness = strandness
    }

    call SamToSortedBam {
        input:
            unsorted_sam = FastqToSam.initially_mapped_sam,
            sample_name = base_file_name
    }

    call PicardMarkDuplicates {
        input:
            sample_name = base_file_name,
            sorted_bam = SamToSortedBam.sorted_bam
    }

    call RemoveUnmappedReads {
        input:
            sample_name = base_file_name,
            duplicates_removed_bam = PicardMarkDuplicates.duplicates_removed_bam,
            single_end_attribute = single_end_attribute,
            remove_unmapped_inline_single_end_script = remove_unmapped_inline_single_end_script,
            remove_unmapped_inline_paired_end_script = remove_unmapped_inline_paired_end_script
    }

    call CallHtseqCount {
        input:
            sample_name = base_file_name,
            sam_to_count = RemoveUnmappedReads.unmapped_reads_removed_sam,
            gtf_annotation = ref_gtf,
            strandness = strandness
    }

    call CompressResults {
        input:
            sample_name = base_file_name,
            htseq_count_txt = CallHtseqCount.htseq_count_txt
    }

    call CollectResultFiles {
        input:
            sorted_bam = SamToSortedBam.sorted_bam,
            htseq_count_compressed_file = CompressResults.htseq_count_compressed_file,
            sorted_bam_result_directory = sorted_bam_result_directory,
            htseq_count_result_directory = htseq_count_result_directory
    }

    # Output files of the workflows.
    output {
        File sorted_bam = SamToSortedBam.sorted_bam
        File htseq_count_txt = CallHtseqCount.htseq_count_txt
    }
}

task StringToFile {
    String sample_name
    String fastq_r1_string
    String fastq_r2_string

    command {
        echo "${fastq_r1_string}" | tr ";" "\n" > ${sample_name}.fastq_r1_list.txt
        echo "${fastq_r2_string}" | tr ";" "\n" > ${sample_name}.fastq_r2_list.txt
    }

    output {
        Array[File] fastq_list_files = glob("*_list.txt")
    }

    runtime {
        memory: "8G"
        cpu: 1
        disks: "local-disk 500 SSD"
        docker: "debian"
    }
}

task FastqToSam {

    Array[File] fastq_list_files

    Array[File]? fastq_r1_list = read_lines(fastq_list_files[0])
    Array[File]? fastq_r2_list = read_lines(fastq_list_files[1])

    Boolean single
    String sample_name

    Array[File]+ hisat_index
    String hisat_prefix

    String? strandness
    String strandness_arg = if defined(strandness) then "--rna-strandness " + strandness + " " else ""

    command {
        if [ "$single" = true ] ; then
            files=$(echo "-U "${sep="," fastq_r1_list})
        else
            files=$(echo "-1 "${sep="," fastq_r1_list}" -2 "${sep="," fastq_r2_list})
        fi

        /usr/local/bin/hisat2 -p 2 --dta -x ${hisat_prefix} ${strandness_arg} $files -S ${sample_name}.initially_mapped.sam
    }

    output {
        File initially_mapped_sam = glob("*.initially_mapped.sam")[0]
    }

    runtime {
        memory: "13G"
        cpu: 2
        disks: "local-disk 500 SSD"
        docker: "zlskidmore/hisat2:latest"
    }
}

task SamToSortedBam {
    File unsorted_sam
    String sample_name

    command {
        /usr/local/bin/samtools sort -@ 2 -l 9 -o ${sample_name}.sorted.bam ${unsorted_sam}
    }

    output {
        File sorted_bam = "${sample_name}.sorted.bam"
    }

    runtime {
        memory: "13G"
        cpu: 2
        disks: "local-disk 500 SSD"
        docker: "zlskidmore/samtools:latest"
    }
}

task PicardMarkDuplicates {
    String sample_name
    File sorted_bam

    command {
        java -jar /usr/picard/picard.jar MarkDuplicates I=${sorted_bam} O=${sample_name}.duplicates_removed.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=${sample_name}.duplicates_removed.txt QUIET=true COMPRESSION_LEVEL=9 VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
    }

    output {
        File duplicates_removed_bam = "${sample_name}.duplicates_removed.bam"
    }

    runtime {
        docker: "broadinstitute/picard:latest"
        disks: "local-disk 500 SSD"
        memory: "13G"
        cpu: 2
    }
}

task RemoveUnmappedReads {
    String sample_name
    File duplicates_removed_bam

    Boolean single_end_attribute
    File remove_unmapped_inline_single_end_script
    File remove_unmapped_inline_paired_end_script

    # whether the sequencing is single-ended or paried-ended determines which Perl script is used to remove unmapped reads from the sam files.
    File remove_unmapped_reads_script = if single_end_attribute then "${remove_unmapped_inline_single_end_script}" else "${remove_unmapped_inline_paired_end_script}"

    command {
        /usr/local/bin/samtools view ${duplicates_removed_bam} | perl ${remove_unmapped_reads_script} > ${sample_name}.unmapped_reads_removed.sam
    }

    output {
        File unmapped_reads_removed_sam = "${sample_name}.unmapped_reads_removed.sam"
    }

    runtime {
        memory: "4G"
        cpu: 1
        disks: "local-disk 500 SSD"
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    }
}

task CallHtseqCount {

    String sample_name

    File sam_to_count
    File gtf_annotation

    # The HTSeq-count needs to know whether the reads are stranded or stranded.
    String? strandness
    String strandness_arg = if defined(strandness) then "--stranded=yes" else "--stranded=no"

    command {
        # htseq-count [options] <alignment_files> <gff_file>
        # https://htseq.readthedocs.io/en/release_0.11.1/count.html
        /usr/local/bin/htseq-count ${strandness_arg} ${sam_to_count} ${gtf_annotation} > ${sample_name}.htseq_count.txt
    }

    output {
        File htseq_count_txt = "${sample_name}.htseq_count.txt"
    }

    runtime {
        # docker: "biocontainers/htseq:v0.11.2-1-deb-py3_cv1"
        docker: "quay.io/biocontainers/htseq:0.11.2--py36h7eb728f_0"
        memory: "8G"
        cpu: 1
        disks: "local-disk 500 SSD"
    }
}

task CompressResults {
    String sample_name
    File htseq_count_txt

    command {
        # gzip -9 -cvf ${htseq_count_txt} > ${sample_name}.htseq_count.txt.gz
        bzip2 -9 -kvf ${htseq_count_txt}
    }

    output {
        # File htseq_count_compressed_file = "${sample_name}.htseq_count.txt.gz"
        File htseq_count_compressed_file = "${sample_name}.htseq_count.txt.bz2"
    }

    runtime {
        memory: "8G"
        cpu: 1
        disks: "local-disk 500 SSD"
        docker: "debian"
    }
}


task CollectResultFiles {
    File sorted_bam
    File htseq_count_compressed_file

    String sorted_bam_result_directory
    String htseq_count_result_directory

    command {
        gsutil cp ${sorted_bam} ${sorted_bam_result_directory}
        gsutil cp ${htseq_count_compressed_file} ${htseq_count_result_directory}
    }

    runtime {
        docker: "google/cloud-sdk:latest"
        memory: "8G"
        cpu: 1
        disks: "local-disk 500 SSD"
    }

}
