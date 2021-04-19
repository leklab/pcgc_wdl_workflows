# This is a sample script for GATK3 workflow wdl script. 
# inputSamplesFile should have 3 columns: sample_id, fastq1, fastq2.
# reference files and gatk/picard are defined as strings to avoid dealing with indices and copying files.
# Make sure that you have bwa, gatk, picard and R ready to be accessed from command line. 
# For gatk and picard input string is specified in inputs.json.
# Tested with the following versions on Ruddle:
# module load BWA/0.7.15-foss-2016a
# module load picard/2.9.0-Java-1.8.0_121
# module load GATK/3.8-0-Java-1.8.0_121
# module load R/3.4.1-foss-2016b


workflow exome_processing {
  
	String ref
	String dbsnp
	String targets

	File inputSamplesFile
	Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)

	String gatk
	String picard

	scatter (sample in inputSamples) {

		call bwa {
			input:
				ref = ref,
				fastq1 = sample[1],
				fastq2 = sample[2],
				prefix = sample[0]
		}

		call sort_sam {
			input:
				picard = picard,
				prefix = sample[0],
				input_sam = bwa.output_sam
		}

		call markdups{
			input:
				picard = picard,
				prefix = sample[0],
				input_bam = sort_sam.output_bam,
		}
	  
		call baserecal_pre {
			input:
				gatk = gatk,
				ref = ref,
				prefix = sample[0],
				input_bam = markdups.output_bam,
				input_bam_index = markdups.output_bam_index,
				dbsnp = dbsnp
		}

		call baserecal_post {
			input:
				gatk = gatk,
				ref = ref,
				prefix = sample[0],
				dbsnp = dbsnp,
				input_bam = markdups.output_bam,
				input_bam_index = markdups.output_bam_index,
				recal_data = baserecal_pre.recal_data

		}
	  
		call anal_covars {
			input:
				gatk = gatk,
				ref = ref,
				prefix = sample[0],
				recal_data = baserecal_pre.recal_data,
				recal_data_post = baserecal_post.recal_data_post

		}

		call print_reads {
			input:
				gatk = gatk,
				ref = ref,
				prefix = sample[0],
				recal_data = baserecal_pre.recal_data,
				input_bam = markdups.output_bam,
				input_bam_index = markdups.output_bam_index

		}

		call haplotypecaller {
			input:
				gatk = gatk,
				ref = ref,
				prefix = sample[0],
				input_bam = print_reads.output_bam,
				input_bam_index = print_reads.output_bam_index,
				dbsnp = dbsnp,
				targets = targets,	  		
		} 

		output {
			Array[File] gvcf = haplotypecaller.output_gvcf
			Array[File] gvcf_index = haplotypecaller.output_gvcf_index
			Array[File] bam = print_reads.output_bam
			Array[File] bam_index = print_reads.output_bam_index
			Array[File] dedup_metrics = markdups.metrics
			Array[File] recal_plot = anal_covars.recal_plots
		}
	}
}

task bwa{
	String fastq1
	String fastq2
	String ref
	String prefix

	command {
		bwa mem -M -t 12 -R "@RG\tID:${prefix}\tPL:ILLUMINA\tSM:${prefix}" ${ref} ${fastq1} ${fastq2} > ${prefix}.aligned_reads.sam 2> ${prefix}.bwa.stderr.log
	}

	runtime {
		cpus: 12
		requested_memory: 16000 
	}

	output {
		File output_sam = "${prefix}.aligned_reads.sam"
		File bwa_stderr_log = "${prefix}.bwa.stderr.log"
	}
}

task sort_sam{
	String picard
	String prefix
	File input_sam

	command {
		java -Xmx14G -jar ${picard} SortSam INPUT=${input_sam} OUTPUT=${prefix}.sorted_reads.bam SORT_ORDER=coordinate
	}

	runtime {
		cpus: 8
		requested_memory: 16000 
	}

	output {
		File output_bam = "${prefix}.sorted_reads.bam"
	}
}

task markdups{
	String picard
	String prefix
	File input_bam

	command {
		java -Xmx14G -jar ${picard} MarkDuplicates INPUT=${input_bam} OUTPUT=${prefix}.dedup_reads.bam METRICS_FILE=${prefix}.metrics.txt CREATE_INDEX=true
	}

	runtime {
		cpus: 8
		requested_memory: 16000 
	}

	output {
		File output_bam = "${prefix}.dedup_reads.bam"
		File metrics = "${prefix}.metrics.txt"
		File output_bam_index = "${prefix}.dedup_reads.bai"
	}
}

task baserecal_pre{
	String gatk
	String ref
	String prefix
	String dbsnp
	File input_bam
	File input_bam_index

	command {
		java -Xmx14G -jar ${gatk} -T BaseRecalibrator -R ${ref} -I ${input_bam} -knownSites ${dbsnp} -o ${prefix}.recal_data.table
	}

	runtime {
		cpus: 8
		requested_memory: 16000 
	}

	output {
		File recal_data = "${prefix}.recal_data.table"
	}

}

task baserecal_post{
	String gatk
	String ref
	String prefix
	String dbsnp
	File input_bam
	File input_bam_index
	File recal_data

	command {
		java -Xmx14G -jar ${gatk} -T BaseRecalibrator -R ${ref} -I ${input_bam} -knownSites ${dbsnp} -BQSR ${recal_data} -o ${prefix}.post_recal_data.table
	}

	runtime {
		cpus: 8
		requested_memory: 16000 
	}

	output {
		File recal_data_post = "${prefix}.post_recal_data.table"
	}
}

task anal_covars{
	String gatk
	String ref
	String prefix
	File recal_data
	File recal_data_post


	command {
		java -Xmx14G -jar ${gatk} -T AnalyzeCovariates -R ${ref} -before ${recal_data} -after ${recal_data_post} -plots ${prefix}.recalibration_plots.pdf
	}

	runtime {
		cpus: 8
		requested_memory: 16000 
	}

	output {
		File recal_plots = "${prefix}.recalibration_plots.pdf"
	}
}


task print_reads{
	String gatk
	String ref
	String prefix
	File recal_data
	File input_bam
	File input_bam_index


	command {
		java -Xmx14G -jar ${gatk} -T PrintReads -R ${ref} -I ${input_bam} -BQSR ${recal_data} -o ${prefix}.bam
	}

	runtime {
		cpus: 8
		requested_memory: 16000 
	}

	output {
		File output_bam = "${prefix}.bam"
		File output_bam_index = "${prefix}.bai"
	}
}

task haplotypecaller{
	String gatk
	String ref
	String prefix
	String dbsnp
	String targets
	File input_bam
	File input_bam_index


	command {
		java -Xmx14G -jar ${gatk} -T HaplotypeCaller -R ${ref} -I ${input_bam} --emitRefConfidence GVCF -L ${targets} -ip 50 --dbsnp ${dbsnp} -o ${prefix}.g.vcf.gz
	}

	runtime {
		cpus: 8
		requested_memory: 16000 
	}

	output {
		File output_gvcf = "${prefix}.g.vcf.gz"
		File output_gvcf_index = "${prefix}.g.vcf.gz.tbi"
	}
}
