nextflow.enable.dsl=2

process TRIM {
  cpus 1

  input:
    path R1
    path R2
    path adapters_file

  output:
    path "trimmed_R1.fastq.gz", emit: R1
    path "trimmed_R2.fastq.gz", emit: R2
  
  script:
    """
    echo \${TRIMMOMATIC}
    if [ `cat ${adapters_file} | wc -l` -gt 0 ]
    then
      java -jar \${TRIMMOMATIC} PE -phred33 -threads ${task.cpus} ${R1} ${R2} trimmed_R1.fastq.gz trimmed_unpaired_R1.fastq.gz trimmed_R2.fastq.gz trimmed_unpaired_R2.fastq.gz ILLUMINACLIP:${adapters_file}:2:40:14:3:true TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36
    else
      echo "Skipping trimming"
      cp ${R1} trimmed_R1.fastq.gz
      cp ${R2} trimmed_R2.fastq.gz
    fi
    """
}

process DEMULTIPLEX {
    input:
      val  sample_list
      path R1
      path R2
    
    output:
      path 'sample_*/*R1.fastq.gz', emit: R1
      path 'sample_*/*R2.fastq.gz', emit: R2

    errorStrategy 'ignore'

    script:
      """
      solo_options=\"\"
      if ${params.solo}
      then
        solo_options=\"--solo --barcode-spacer ${params.barcode_spacer}\"
      fi
      rust-demultiplex \${solo_options} --umi-group-length ${params.umi_group_length} --umi-length ${params.umi_length} --pcr-primer ${params.pcr_primer} <(gunzip -c ${R1}) <(gunzip -c ${R2}) /dev/null ${sample_list}
      """
}

process COLLAPSE_AND_ASSEMBLE {
    cpus 1

    input:
      tuple val(sample), val(umi_bin), path(R1), path(R2)
    
    output:
      tuple val(sample), val(umi_bin), path("${sample}_${umi_bin}_output.fa.gz"), emit: fasta
      tuple val(sample), val(umi_bin), path("${sample}_${umi_bin}_output.csv.gz"), emit: csv

    script:
      """
      rust_solo_options=\"\"
      assembly_solo_options=\"\"
      if ${params.solo}
      then
        rust_solo_options=\" --solo --barcode-spacer ${params.barcode_spacer} \"
        assembly_solo_options=\" --solo-contig-barcodes ${params.solo_contig_barcodes} --solo-barcodes ${params.solo_barcodes} --assembly-iterations ${params.assembly_iterations} --sampling-target ${params.sampling_target} \"
      fi

      rust-umi \${rust_solo_options} --min-count ${params.min_count} --umi-length ${params.umi_length} --pcr-primer ${params.pcr_primer} <(gunzip -c ${R1}) <(gunzip -c ${R2}) /dev/null ${sample}
      spades_assembler.py \${assembly_solo_options} --r1-orientation ${params.r1_orientation} --threads ${task.cpus} --pcr-primer ${params.pcr_primer} --anchor-start ${params.anchor_start} --anchor-end ${params.anchor_end} output ${sample}_${umi_bin}
      rm -r -f output
      for file in *.fa *.csv
      do
        pigz -p ${task.cpus} \${file}
      done
      """
}

 process COLLECT_RESULTS {
    input:
      tuple val(sample), path(fasta), path(csv)

    output:
      tuple val(sample), path("output/${sample}.fa.gz"), emit: fasta
      tuple val(sample), path("output/${sample}.csv.gz"), emit: csv

    publishDir "${params.out_dir}", mode: 'copy', overwrite: false

    script:
        """
        mkdir -p output
        cat *.fa.gz > output/${sample}.fa.gz
        cat *.csv.gz > output/${sample}.csv.gz
        """
 }

workflow {
    r1 = Channel.from(file(params.R1))
    r2 = Channel.from(file(params.R2))
    adapters_file = Channel.from(file(params.adapters_file))
    
    sample_list = Channel.from(params.sample_list).flatMap {
      arg -> arg.split(',')
    }

    trim_result = TRIM(r1, r2, adapters_file)

    demultiplex_result = DEMULTIPLEX(sample_list, trim_result.R1.first(), trim_result.R2.first())

    assembly_result = COLLAPSE_AND_ASSEMBLE(demultiplex_result.R1.flatten().merge(demultiplex_result.R2.flatten()) {
      r1, r2 -> 
        splat = r1.toString().split("/")[-1].split("_")
        sample = splat[0]
        umi_bin = splat[1]
        return [sample, umi_bin, r1, r2] 
    })
    
    COLLECT_RESULTS(assembly_result.fasta
      .join(assembly_result.csv, by: [0, 1])
      .groupTuple()
      .map {
        sample, umi_bin, fasta, csv -> 
        return [sample, fasta, csv]
      })
    
}