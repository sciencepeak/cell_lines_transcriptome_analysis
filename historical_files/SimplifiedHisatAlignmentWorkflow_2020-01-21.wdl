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
    Boolean? single_ended
    Boolean single_default = select_first([single_ended, false])

    # Process input args: usually the library prepation for strandness is "RF".
    String? strandness
    String strandness_arg = if defined(strandness) then "--rna-strandness " + strandness + " " else ""

    # REFERENCE FILES for HiSAT2
    File hisat_index_file
    Array[String] hisat_index = read_lines(hisat_index_file)
    String hisat_prefix

    # Reference Files for BaseRecalibrator and/or ApplyBQSR
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # GATK4 Reference Files for BaseRecalibrator and/or ApplyBQSR
    File known_indels
    File known_indels_index
    File thousG
    File thousG_index
    File dbsnp
    File dbsnp_index


    call FastqToSam {
        input:
            fastq_1 = fastqr1,
            fastq_2 = fastqr2,
            single = single_default,
            sample_name = base_file_name,
            hisat_index = hisat_index,
            hisat_prefix = hisat_prefix,
            strandness = strandness_arg
    }

    call SamToSortedBam {
        input:
            sam = FastqToSam.mapped_sam,
            sample_name = base_file_name,
            id = FastqToSam.id,
            sm = FastqToSam.sm
    }

    call PicardMarkDuplicates {
        input:
            sample_name = base_file_name,
            id = FastqToSam.id,
            sm = FastqToSam.sm,
            sorted_bam_file = SamToSortedBam.sorted_bam
    }

    call runPicard_Gatk {
        input:
            sample_name = base_file_name,
            id = FastqToSam.id,
            sm = FastqToSam.sm,
            bamFile = PicardMarkDuplicates.md_bam,

            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,

            known_indels =known_indels,
            known_indels_index = known_indels_index,
            thousG = thousG,
            thousG_index = thousG_index,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index
    }


    # Output BAM files
    output {
        File mapped_sam = FastqToSam.mapped_sam
        File sorted_bam = SamToSortedBam.sorted_bam
        File md_bam = PicardMarkDuplicates.md_bam
        File finalBam = runPicard_Gatk.finalBam
    }
}

task FastqToSam {
    File fastq_1
    File? fastq_2
    Boolean single
    String sample_name

    Array[File]+ hisat_index
    String hisat_prefix
    String strandness

    command <<<
        if [ "$single" = true ] ; then
            files=$(echo "-U "${fastq_1})
        else
            files=$(echo "-1 "${fastq_1}" -2 "${fastq_2})
        fi

        id=$(zcat < ${fastq_1} | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/./g')
        sm=$(zcat < ${fastq_1} | head -n 1 | grep -Eo "[ATGCN]+$")

        echo $id > ID.txt
        echo $sm > SM.txt

        /usr/local/bin/hisat2 -p 8 --dta -x ${hisat_prefix} --rg-id $id --rg PL:ILLUMINA --rg PU:${sample_name} --rg LB:$id.$sm --rg SM:${sample_name} ${strandness} $files -S ${sample_name}.$id.$sm.aligned.sam
    >>>

    output {
        File mapped_sam = glob("*.aligned.sam")[0]
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

    command <<<
        /usr/local/bin/samtools sort -@ -8 -o ${sample_name}.${id}.${sm}.sorted.bam ${sam}
    >>>

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
    String sm
    String id
    File sorted_bam_file

    command <<<
        java -Xmx4g -jar /usr/picard/picard.jar MarkDuplicates I=${sorted_bam_file} O=${sample_name}.${id}.${sm}.md.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=${sample_name}.${id}.${sm}.md.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
    >>>

    output {
        File md_bam = "${sample_name}.${id}.${sm}.md.bam"
    }

    runtime {
        docker: "broadinstitute/picard:latest"
        disks: "local-disk 500 SSD"
        memory: "16G"
        cpu: 1
    }
}

task runPicard_Gatk {
    String sample_name
    String sm
    String id
    File bamFile

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    File known_indels
    File known_indels_index
    File thousG
    File thousG_index
    File dbsnp
    File dbsnp_index

    command <<<
        /usr/gitc/gatk4/gatk-launch BaseRecalibrator -R ${ref_fasta} -I ${bamFile} -O ${sample_name}.bqsr.table -knownSites ${thousG} -knownSites ${known_indels} -knownSites ${dbsnp}

        /usr/gitc/gatk4/gatk-launch ApplyBQSR -R ${ref_fasta} -I ${bamFile} -O ${sample_name}.FINAL.bam -bqsr ${sample_name}.bqsr.table --static_quantized_quals 10 --static_quantized_quals 20 --static_quantized_quals 30
    >>>

    output {
        File finalBam = "${sample_name}.FINAL.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 500 SSD"
        memory: "16G"
        cpu: 1
    }
}
