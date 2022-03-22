process art_illumina {

    tag { assembly_id + ' / ' + fold_coverage + 'x' + ' / ' + 'len=' + read_length + ' / ' + 'replicate=' + replicate }

    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}*_read_simulation_parameters.csv", mode: 'copy'

    input:
    tuple val(assembly_id), path(assembly), val(fold_coverage), val(replicate)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}*_R1.fastq"), path("${assembly_id}-${md5_fragment}*_R2.fastq"), emit: reads
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_read_simulation_parameters.csv"), emit: metrics

    script:
    mean_fragment_length = params.mean_fragment_length
    stdev_fragment_length = params.stdev_fragment_length
    quality_shift_r1 = params.quality_shift_r1
    quality_shift_r2 = params.quality_shift_r2
    read_length = params.read_length
    seed = Math.round(Math.random() * 1000000)
    md5_input = assembly_id + fold_coverage.toString() + read_length.toString() + mean_fragment_length.toString() + stdev_fragment_length.toString() + seed.toString() + quality_shift_r1.toString() + quality_shift_r2.toString()
    md5_fragment = md5_input.md5()[0..3]
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    art_illumina \
      --paired \
      --in ${assembly} \
      --fcov ${fold_coverage} \
      --len ${read_length} \
      --mflen ${mean_fragment_length} \
      --sdev ${stdev_fragment_length} \
      --rndSeed ${seed} \
      --qShift ${params.quality_shift_r1} \
      --qShift2 ${params.quality_shift_r2} \
      --out ${assembly_id}-${md5_fragment}_R
    mv ${assembly_id}-${md5_fragment}_R1.fq ${assembly_id}-${md5_fragment}_R1.fastq
    mv ${assembly_id}-${md5_fragment}_R2.fq ${assembly_id}-${md5_fragment}_R2.fastq
    echo 'sample_id,replicate,random_seed,fold_coverage,read_length,mean_fragment_length,stdev_fragment_length,quality_shift_r1,quality_shift_r2' > ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    echo '${assembly_id}-${md5_fragment},${replicate},${seed},${fold_coverage},${read_length},${mean_fragment_length},${stdev_fragment_length},${quality_shift_r1},${quality_shift_r2}' >> ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    """
}

process simulate_contaminant_reads {

    tag { contaminant_id + ' / ' + proportion}

    input:
    tuple val(contaminant_id), path(assembly), val(proportion)

    output:
    tuple val(contaminant_id), path("${contaminant_id}*_R1.fastq"), path("${contaminant_id}*_R2.fastq"), val(proportion)

    script:
    mean_fragment_length = params.mean_fragment_length
    stdev_fragment_length = params.stdev_fragment_length
    quality_shift_r1 = params.quality_shift_r1
    quality_shift_r2 = params.quality_shift_r2
    read_length = params.read_length
    seed = Math.round(Math.random() * 1000000)
    """
    art_illumina \
      --paired \
      --in ${assembly} \
      --fcov 30 \
      --len ${read_length} \
      --mflen ${mean_fragment_length} \
      --sdev ${stdev_fragment_length} \
      --rndSeed ${seed} \
      --qShift ${params.quality_shift_r1} \
      --qShift2 ${params.quality_shift_r2} \
      --out ${contaminant_id}_R
    mv ${contaminant_id}_R1.fq ${contaminant_id}_R1.fastq
    mv ${contaminant_id}_R2.fq ${contaminant_id}_R2.fastq
    """
}

process downsample_simulated_reads {
  tag { assembly_id + '-' + md5_fragment + ' / ' + proportion_uncontaminated }

  input:
  tuple val(assembly_id), val(md5_fragment), path(assembly_reads_r1), path(assembly_reads_r2), val(proportion_contaminants)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}*_R1.fastq"), path("${assembly_id}-${md5_fragment}*_R2.fastq")

  script:
  proportion_uncontaminated = 1.0 - proportion_contaminants
  seed = Math.round(Math.random() * 1000000)
  """
  seqkit sample -s ${seed} -p ${proportion_uncontaminated} ${assembly_reads_r1} -o ${assembly_id}-${md5_fragment}_sample_R1.fastq
  seqkit sample -s ${seed} -p ${proportion_uncontaminated} ${assembly_reads_r2} -o ${assembly_id}-${md5_fragment}_sample_R2.fastq
  """
}

process downsample_contaminant_reads {

  tag { assembly_id + '-' + md5_fragment + ' / ' + contaminant_id + ' / ' + contaminant_proportion }

  publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${contaminant_id}_sample_R*.fastq.gz", mode: 'copy'
  publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv", mode: 'copy'

  input:
  tuple val(contaminant_id), path(contaminant_reads_r1), path(contaminant_reads_r2), val(contaminant_proportion), val(assembly_id), val(md5_fragment), path(assembly_reads_r1), path(assembly_reads_r2)

  output:
  tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}_sample_R1.fastq"), path("${contaminant_id}_sample_R2.fastq"), emit: uncompressed_reads
  tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}_sample_R1.fastq.gz"), path("${contaminant_id}_sample_R2.fastq.gz"), emit: compressed_reads
  tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv"), emit: num_reads_csv

  script:
  output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
  seed = Math.round(Math.random() * 1000000)
  """
  seqkit stats -T ${assembly_reads_r1} | tail -n 1 | cut -f 4 | tr -d ',' > num_simulated_read_pairs
  echo 'sample_id,contaminant_id,num_simulated_read_pairs,num_contaminant_read_pairs,target_contaminant_proportion' > ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv
  python -c "import sys; print(int(round(int(sys.stdin.read().strip()) * ${contaminant_proportion})))" < num_simulated_read_pairs > num_contaminant_read_pairs
  paste -d ',' <(echo "${assembly_id}-${md5_fragment}") <(echo "${contaminant_id}") num_simulated_read_pairs num_contaminant_read_pairs <(echo "${contaminant_proportion}") >> ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv
  seqkit sample -s ${seed} -n \$(cat num_contaminant_read_pairs) ${contaminant_reads_r1} > ${contaminant_id}_sample_R1.fastq
  seqkit sample -s ${seed} -n \$(cat num_contaminant_read_pairs) ${contaminant_reads_r2} > ${contaminant_id}_sample_R2.fastq
  gzip --keep ${contaminant_id}_sample_R*.fastq
  """
}

process introduce_contaminants {

  tag { assembly_id + '-' + md5_fragment }

  input:
  tuple val(assembly_id), val(md5_fragment), path(assembly_r1), path(assembly_r2), val(contaminant_ids), path(contaminants_r1), path(contaminants_r2)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_R1.fastq"), path("${assembly_id}-${md5_fragment}_R2.fastq")

  script:
  seed = Math.round(Math.random() * 1000000)
  """
  mv ${assembly_r1} uncontaminated_R1.fastq
  mv ${assembly_r2} uncontaminated_R2.fastq
  cat uncontaminated_R1.fastq ${contaminants_r1} > ${assembly_id}-${md5_fragment}_unshuffled_R1.fastq
  cat uncontaminated_R2.fastq ${contaminants_r2} > ${assembly_id}-${md5_fragment}_unshuffled_R2.fastq
  paste <(cat ${assembly_id}-${md5_fragment}_unshuffled_R1.fastq) <(cat ${assembly_id}-${md5_fragment}_unshuffled_R2.fastq) \
    | paste - - - - | shuf | awk -F'\\t' '{OFS="\\n"; print \$1,\$3,\$5,\$7 > "${assembly_id}-${md5_fragment}_R1.fastq"; print \$2,\$4,\$6,\$8 > "${assembly_id}-${md5_fragment}_R2.fastq"}'
  """
}

process fastp {

  tag { assembly_id + '-' + md5_fragment }

  publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_R{1,2}.fastq.gz", mode: 'copy'
  publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_fastp.json", mode: 'copy'

  input:
  tuple val(assembly_id), val(md5_fragment), path(reads_1), path(reads_2)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.json"), emit: json
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.csv"), emit: csv
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_R{1,2}.fastq.gz"), emit: untrimmed_reads

  script:
  output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
  """
  fastp -i ${reads_1} -I ${reads_2} -o ${assembly_id}-${md5_fragment}_trimmed_R1.fastq.gz -O ${assembly_id}-${md5_fragment}_trimmed_R2.fastq.gz
  mv fastp.json ${assembly_id}-${md5_fragment}_fastp.json
  fastp_json_to_csv.py -s ${assembly_id}-${md5_fragment} ${assembly_id}-${md5_fragment}_fastp.json > ${assembly_id}-${md5_fragment}_fastp.csv
  cp ${reads_1} untrimmed_R1.fastq
  cp ${reads_2} untrimmed_R2.fastq
  gzip -c untrimmed_R1.fastq > ${assembly_id}-${md5_fragment}_R1.fastq.gz
  gzip -c untrimmed_R2.fastq > ${assembly_id}-${md5_fragment}_R2.fastq.gz
  """
}

process bwa_align {

  tag { assembly_id + '-' + md5_fragment + ' / ' + assembly_id }

  publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}.{bam,bam.bai}", mode: 'copy'

  input:
  tuple val(assembly_id), val(md5_fragment), path(reads_1), path(reads_2), path(ref)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}.bam"), path("${assembly_id}-${md5_fragment}.bam.bai")

  script:
  output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
  """
  bwa index ${ref}
  bwa mem -t ${task.cpus} ${ref} ${reads_1} ${reads_2} | \
    samtools sort -o ${assembly_id}-${md5_fragment}.bam -
  samtools index ${assembly_id}-${md5_fragment}.bam
  """
}

process qualimap_bamqc {

  tag { assembly_id + '-' + md5_fragment }

  publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_bamqc"

  input:
  tuple val(assembly_id), val(md5_fragment), file(alignment), file(alignment_index)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_bamqc/genome_results.txt"), emit: genome_results
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_bamqc"), emit: bamqc_dir
  
  script:
  output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
  """
  qualimap bamqc -bam ${alignment} --outdir ${assembly_id}-${md5_fragment}_bamqc
  """
}

process qualimap_bamqc_genome_results_to_csv {

  tag { assembly_id + '-' + md5_fragment }

  publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_qualimap_bamqc_genome_results.csv"

  executor 'local'

  input:
  tuple val(assembly_id), val(md5_fragment), path(qualimap_bamqc_genome_results)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_qualimap_bamqc_genome_results.csv")

  script:
  output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
  """
  qualimap_bamqc_genome_results_to_csv.py -s ${assembly_id}-${md5_fragment} ${qualimap_bamqc_genome_results} > ${assembly_id}-${md5_fragment}_qualimap_bamqc_genome_results.csv
  """
}

process samtools_stats {

  tag { assembly_id + '-' + md5_fragment }

  publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"

  input:
  tuple val(assembly_id), val(md5_fragment), file(alignment), file(alignment_index)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"), emit: summary
  
  script:
  output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
  """
  samtools stats ${alignment} > ${assembly_id}-${md5_fragment}_samtools_stats.txt
  grep ^SN ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2- > ${assembly_id}-${md5_fragment}_samtools_stats_summary.txt
  """
}

process mosdepth {

  tag { assembly_id + '-' + md5_fragment }

  publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"

  input:
  tuple val(assembly_id), val(md5_fragment), file(alignment), file(alignment_index), val(depth_by)

  output:
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_by_${depth_by}.mosdepth.summary.txt"), emit: summary
  tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_by_${depth_by}.mosdepth.summary.txt"), emit: regions
  
  script:
  output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
  """
  mosdepth -t ${task.cpus} --fast-mode --by ${depth_by} --no-per-base ${assembly_id}-${md5_fragment}_by_${depth_by} ${alignment}
  gunzip ${assembly_id}-${md5_fragment}_by_${depth_by}.regions.bed.gz
  mv ${assembly_id}-${md5_fragment}_by_${depth_by}.regions.bed ${assembly_id}-${md5_fragment}_by_${depth_by}.mosdepth.regions.bed
  """
}