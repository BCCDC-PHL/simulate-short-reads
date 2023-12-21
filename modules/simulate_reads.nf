process art_illumina {

    tag { assembly_id + ' / ' + fold_coverage + 'x' + ' / ' + 'len=' + read_length + ' / ' + 'replicate=' + replicate }

    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}*_read_simulation_parameters.csv", mode: 'copy'

    input:
    tuple val(assembly_id), path(assembly), val(fold_coverage), val(replicate)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}_R1.fastq.gz"), path("${assembly_id}_R2.fastq.gz"), emit: reads
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
      --out ${assembly_id}_R

    mv ${assembly_id}_R1.fq ${assembly_id}_R1.fastq
    mv ${assembly_id}_R2.fq ${assembly_id}_R2.fastq

    gzip ${assembly_id}_R1.fastq
    gzip ${assembly_id}_R2.fastq

    rm -f ${assembly_id}-${md5_fragment}_R1.aln
    rm -f ${assembly_id}-${md5_fragment}_R2.aln
    
    echo 'sample_id,replicate,random_seed,fold_coverage,read_length,mean_fragment_length,stdev_fragment_length,quality_shift_r1,quality_shift_r2' > ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    echo '${assembly_id}-${md5_fragment},${replicate},${seed},${fold_coverage},${read_length},${mean_fragment_length},${stdev_fragment_length},${quality_shift_r1},${quality_shift_r2}' >> ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    """
}

process simulate_contaminant_reads {

    tag { assembly_id + '-' + md5_fragment + ' / ' + contaminant_id + ' / ' + proportion }

    input:
    tuple val(contaminant_id), path(assembly), val(proportion), val(assembly_id), val(md5_fragment), path(reads_r1), path(reads_r2)

    output:
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}*_R1.fastq.gz"), path("${contaminant_id}*_R2.fastq.gz"), val(proportion)

    script:
    mean_fragment_length = params.mean_fragment_length
    stdev_fragment_length = params.stdev_fragment_length
    quality_shift_r1 = params.quality_shift_r1
    quality_shift_r2 = params.quality_shift_r2
    read_length = params.read_length
    seed = Math.round(Math.random() * 1000000)
    """
    contaminant_genome_size=\$(seqkit stats --tabular ${assembly} | cut -f 5 | tail -n 1)
    num_simulated_reads_r1=\$(seqkit stats --tabular ${reads_r1} | cut -f 4 | tail -n 1)
    num_simulated_reads_r2=\$(seqkit stats --tabular ${reads_r2} | cut -f 4 | tail -n 1)
    num_simulated_reads=\$(( num_simulated_reads_r1 + num_simulated_reads_r2))
    fcov=\$(calculate_required_coverage.py --num-simulated-reads \$num_simulated_reads --read-length ${read_length} --contaminant-genome-size \$contaminant_genome_size --contaminant-proportion ${proportion})

    art_illumina \
      --paired \
      --in ${assembly} \
      --fcov \$fcov \
      --len ${read_length} \
      --mflen ${mean_fragment_length} \
      --sdev ${stdev_fragment_length} \
      --rndSeed ${seed} \
      --qShift ${params.quality_shift_r1} \
      --qShift2 ${params.quality_shift_r2} \
      --out ${contaminant_id}_R

    mv ${contaminant_id}_R1.fq ${contaminant_id}_R1.fastq
    mv ${contaminant_id}_R2.fq ${contaminant_id}_R2.fastq

    gzip ${contaminant_id}_R1.fastq
    gzip ${contaminant_id}_R2.fastq
    """
}

process downsample_simulated_reads {

    tag { assembly_id + '-' + md5_fragment + ' / ' + proportion_uncontaminated }

    input:
    tuple val(assembly_id), val(md5_fragment), path(assembly_reads_r1), path(assembly_reads_r2), val(proportion_contaminants)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-sample_R1.fastq.gz"), path("${assembly_id}-sample_R2.fastq.gz")

    script:
    proportion_uncontaminated = 1.0 - proportion_contaminants
    seed = Math.round(Math.random() * 1000000)
    """
    seqkit sample -s ${seed} -p ${proportion_uncontaminated} ${assembly_reads_r1} -o ${assembly_id}-sample_R1.fastq
    seqkit sample -s ${seed} -p ${proportion_uncontaminated} ${assembly_reads_r2} -o ${assembly_id}-sample_R2.fastq

    gzip ${assembly_id}-sample_R1.fastq
    gzip ${assembly_id}-sample_R2.fastq
    """
}

process downsample_contaminant_reads {

    tag { assembly_id + '-' + md5_fragment + ' / ' + contaminant_id + ' / ' + contaminant_proportion }

    publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${contaminant_id}_contaminant_R*.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv", mode: 'copy'

    input:
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path(contaminant_reads_r1), path(contaminant_reads_r2), val(contaminant_proportion), path(assembly_reads_r1), path(assembly_reads_r2)

    output:
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}_contaminant_R1.fastq.gz"), path("${contaminant_id}_contaminant_R2.fastq.gz"), emit: reads
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv"), emit: num_reads_csv

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    seed = Math.round(Math.random() * 1000000)
    """
    seqkit stats -T ${assembly_reads_r1} | tail -n 1 | cut -f 4 | tr -d ',' > num_simulated_read_pairs
    echo 'sample_id,contaminant_id,num_simulated_read_pairs,num_contaminant_read_pairs,target_contaminant_proportion' > ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv
    python -c "import sys; print(int(round(int(sys.stdin.read().strip()) * ${contaminant_proportion})))" < num_simulated_read_pairs > num_contaminant_read_pairs
    paste -d ',' <(echo "${assembly_id}-${md5_fragment}") <(echo "${contaminant_id}") num_simulated_read_pairs num_contaminant_read_pairs <(echo "${contaminant_proportion}") >> ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv
    seqkit sample -s ${seed} -n \$(cat num_contaminant_read_pairs) ${contaminant_reads_r1} | gzip > ${contaminant_id}_contaminant_R1.fastq.gz
    seqkit sample -s ${seed} -n \$(cat num_contaminant_read_pairs) ${contaminant_reads_r2} | gzip > ${contaminant_id}_contaminant_R2.fastq.gz
    """
}

process introduce_contaminants {

    tag { assembly_id + '-' + md5_fragment + ' / ' + contaminant_ids }

    input:
    tuple val(assembly_id), val(md5_fragment), path(assembly_r1), path(assembly_r2), val(contaminant_ids), path(contaminants_r1), path(contaminants_r2)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-contaminated_R1.fastq"), path("${assembly_id}-contaminated_R2.fastq"), emit: reads

    script:
    seed = Math.round(Math.random() * 1000000)
    """
    cp ${assembly_r1} uncontaminated_R1.fastq.gz
    cp ${assembly_r2} uncontaminated_R2.fastq.gz

    gunzip uncontaminated_R1.fastq.gz
    gunzip uncontaminated_R2.fastq.gz

    cp ${contaminants_r1} contaminants_R1.fastq.gz
    cp ${contaminants_r2} contaminants_R2.fastq.gz

    gunzip contaminants_R1.fastq.gz
    gunzip contaminants_R2.fastq.gz
    
    cat uncontaminated_R1.fastq contaminants_R1.fastq > ${assembly_id}-${md5_fragment}_unshuffled_R1.fastq
    cat uncontaminated_R2.fastq contaminants_R1.fastq > ${assembly_id}-${md5_fragment}_unshuffled_R2.fastq
    paste <(cat ${assembly_id}-${md5_fragment}_unshuffled_R1.fastq) <(cat ${assembly_id}-${md5_fragment}_unshuffled_R2.fastq) \
	| paste - - - - | shuf | awk -F'\\t' '{OFS="\\n"; print \$1,\$3,\$5,\$7 > "${assembly_id}-contaminated_R1.fastq"; print \$2,\$4,\$6,\$8 > "${assembly_id}-contaminated_R2.fastq"}'
    """
}

process fastp {

    tag { assembly_id + '-' + md5_fragment }

    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_R{1,2}.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_fastp.{json,csv}", mode: 'copy'

    input:
    tuple val(assembly_id), val(md5_fragment), path(reads_1), path(reads_2)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.json"), emit: json
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.csv"), emit: fastp_csv
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_R{1,2}.fastq.gz"), emit: untrimmed_reads  

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    fastp \
	-i ${reads_1} \
	-I ${reads_2}

    mv fastp.json ${assembly_id}-${md5_fragment}_fastp.json

    fastp_json_to_csv.py \
	-s ${assembly_id}-${md5_fragment} \
	${assembly_id}-${md5_fragment}_fastp.json \
	> ${assembly_id}-${md5_fragment}_fastp.csv

    mv ${reads_1} ${assembly_id}-${md5_fragment}_R1.fastq.gz
    mv ${reads_2} ${assembly_id}-${md5_fragment}_R2.fastq.gz
    """
}

process bwa_align {

    tag { assembly_id + '-' + md5_fragment + ' / ' + assembly_id }

    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}.{bam,bam.bai}", mode: 'copy', enabled: params.keep_bams

    input:
    tuple val(assembly_id), val(md5_fragment), path(reads_1), path(reads_2), path(ref)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}.bam"), path("${assembly_id}-${md5_fragment}.bam.bai")

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    bwa index ${ref}

    bwa mem \
	-t ${task.cpus} \
	${ref} \
	${reads_1} \
	${reads_2} \
	| samtools sort -o ${assembly_id}-${md5_fragment}.bam
    
    samtools index ${assembly_id}-${md5_fragment}.bam
    """
}

process qualimap_bamqc {

    tag { assembly_id + '-' + md5_fragment }

    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_bamqc"
    publishDir  "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_short_qualimap_alignment_qc.csv"

    input:
    tuple val(assembly_id), val(md5_fragment), file(alignment), file(alignment_index)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_bamqc/genome_results.txt"), emit: genome_results
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_bamqc"), emit: bamqc_dir
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_qualimap_alignment_qc.csv"), emit: alignment_qc
  
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    qualimap bamqc \
	-bam ${alignment} \
	--outdir ${assembly_id}-${md5_fragment}_bamqc

    qualimap_bamqc_genome_results_to_csv.py \
	-s ${assembly_id}-${md5_fragment} \
	${assembly_id}-${md5_fragment}_bamqc/genome_results.txt \
	> ${assembly_id}-${md5_fragment}_qualimap_alignment_qc.csv
    """
}


process samtools_stats {

    tag { assembly_id + '-' + md5_fragment }

    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_samtools_stats*"

    input:
    tuple val(assembly_id), val(md5_fragment), path(alignment), path(alignment_index)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats.txt"), emit: stats
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"), emit: stats_summary
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_summary.csv"), emit: stats_summary_csv
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_insert_sizes.tsv"), emit: insert_sizes
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_coverage_distribution.tsv"), emit: coverage_distribution

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    samtools stats \
	--threads ${task.cpus} \
	${alignment[0]} > ${assembly_id}-${md5_fragment}_samtools_stats.txt

    grep '^SN' ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2-  > ${assembly_id}-${md5_fragment}_samtools_stats_summary.txt

    parse_samtools_stats_summary.py -i ${assembly_id}-${md5_fragment}_samtools_stats_summary.txt -s ${assembly_id} > ${assembly_id}-${md5_fragment}_samtools_stats_summary.csv

    echo "insert_size,pairs_total,inward_oriented_pairs,outward_oriented_pairs,other_pairs" | tr ',' '\t' > ${assembly_id}-${md5_fragment}_samtools_stats_insert_sizes.tsv
    grep '^IS' ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2-  >> ${assembly_id}-${md5_fragment}_samtools_stats_insert_sizes.tsv

    echo "coverage,depth" | tr ',' '\t' > ${assembly_id}-${md5_fragment}_samtools_stats_coverage_distribution.tsv
    grep '^COV' ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2- >> ${assembly_id}-${md5_fragment}_samtools_stats_coverage_distribution.tsv	
    """
}


process combine_alignment_qc {

    tag { assembly_id + '-' + md5_fragment }

    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_combined_alignment_qc.csv"

    input:
    tuple val(assembly_id), val(md5_fragment), path(qualimap_genome_results_csv), path(samtools_stats_summary_csv)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_combined_alignment_qc.csv")

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    combine_alignment_qc.py \
	--sample-id ${assembly_id}-${md5_fragment} \
	--read-type "short" \
	--qualimap-bamqc-genome-results ${qualimap_genome_results_csv} \
	--samtools-stats-summary ${samtools_stats_summary_csv} \
	> ${assembly_id}-${md5_fragment}_combined_alignment_qc.csv
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
