workflow SimplifiedHisatAlignmentWorkflow {
    # IMPORTANT:
    # This script only handles fastq. fastqList is not handled.

    # FILE INFO
    # The base_file_name is actually the sample_name
    String base_file_name

    # Required for fastq
    File? fastqr1
    File? fastqr2

    # Handle single_ended or pair_ended?
    # we can select the first valid element.
    # Given an array of optional values, select_first will select the first defined value and return it.
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
    String picard_remove_duplicates_results
    String htseq_count_results

    call FastqToSam {
        input:
            fastq_1 = fastqr1,
            fastq_2 = fastqr2,
            single = single_end_attribute,
            sample_name = base_file_name,
            hisat_index = hisat_index,
            hisat_prefix = hisat_prefix,
            strandness = strandness
    }

    call SamToSortedBam {
        input:
            sam = FastqToSam.initially_mapped_sam,
            sample_name = base_file_name,
            id = FastqToSam.id,
            sm = FastqToSam.sm
    }

    call PicardMarkDuplicates {
        input:
            sample_name = base_file_name,
            sorted_bam_file = SamToSortedBam.sorted_bam,
            picard_remove_duplicates_results = picard_remove_duplicates_results
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
            strandness = strandness,
            htseq_count_results = htseq_count_results
    }


    call CollectResultFiles {
        input:
            duplicates_removed_bam = PicardMarkDuplicates.duplicates_removed_bam,
            htseq_count_txt = CallHtseqCount.htseq_count_txt,
            picard_remove_duplicates_results = picard_remove_duplicates_results,
            htseq_count_results = htseq_count_results
    }

    # Output files of the workflows.
    output {
        File initially_mapped_sam = FastqToSam.initially_mapped_sam
        File sorted_bam = SamToSortedBam.sorted_bam
        File duplicates_removed_bam = PicardMarkDuplicates.duplicates_removed_bam
        File unmapped_reads_removed_sam = RemoveUnmappedReads.unmapped_reads_removed_sam
        File htseq_count_txt = CallHtseqCount.htseq_count_txt
    }
}

task FastqToSam {
    File fastq_1
    File? fastq_2
    Boolean single
    String sample_name

    Array[File]+ hisat_index
    String hisat_prefix

    String? strandness
    String strandness_arg = if defined(strandness) then "--rna-strandness " + strandness + " " else ""

    # use wdl to de the condition detmination is much more elegant than linux command.
    # Linux commands should be as simple as possible.
    command {
        if [ "$single" = true ] ; then
            files=$(echo "-U "${fastq_1})
        else
            files=$(echo "-1 "${fastq_1}" -2 "${fastq_2})
        fi

        id=$(zcat < ${fastq_1} | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/./g')
        sm=$(zcat < ${fastq_1} | head -n 1 | grep -Eo "[ATGCN]+$")

        echo $id > ID.txt
        echo $sm > SM.txt

        /usr/local/bin/hisat2 -p 8 --dta -x ${hisat_prefix} --rg-id $id --rg PL:ILLUMINA --rg PU:${sample_name} --rg LB:$id.$sm --rg SM:${sample_name} ${strandness_arg} $files -S ${sample_name}.$id.$sm.initially_mapped.sam
    }

    output {
        File initially_mapped_sam = glob("*.initially_mapped.sam")[0]
        String id = read_string("ID.txt")
        String sm = read_string("SM.txt")
    }

    runtime {
        memory: "16G"
        cpu: 1
        disks: "local-disk 500 SSD"
        docker: "zlskidmore/hisat2:latest"
    }
}

task SamToSortedBam {
    File sam
    String sample_name
    String sm
    String id

    command {
        /usr/local/bin/samtools sort -@ 8 -o ${sample_name}.${id}.${sm}.sorted.bam ${sam}
    }

    output {
        File sorted_bam = "${sample_name}.${id}.${sm}.sorted.bam"
    }

    runtime {
        memory: "8G"
        cpu: 1
        disks: "local-disk 500 SSD"
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    }
}

task PicardMarkDuplicates {
    String sample_name
    File sorted_bam_file
    String picard_remove_duplicates_results

    command {
        java -Xmx4g -jar /usr/picard/picard.jar MarkDuplicates I=${sorted_bam_file} O=${sample_name}.duplicates_removed.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=${sample_name}.duplicates_removed.txt QUIET=true COMPRESSION_LEVEL=9 VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
    }

    output {
        File duplicates_removed_bam = "${sample_name}.duplicates_removed.bam"
    }

    runtime {
        docker: "broadinstitute/picard:latest"
        disks: "local-disk 500 SSD"
        memory: "16G"
        cpu: 1
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

    String htseq_count_results

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


task CollectResultFiles {
    File duplicates_removed_bam
    File htseq_count_txt
    String picard_remove_duplicates_results
    String htseq_count_results

    command {
        gsutil cp ${duplicates_removed_bam} ${picard_remove_duplicates_results}
        gsutil cp ${htseq_count_txt} ${htseq_count_results}
    }

    runtime {
        docker: "google/cloud-sdk:latest"
        memory: "8G"
        cpu: 1
        disks: "local-disk 500 SSD"
    }

}
