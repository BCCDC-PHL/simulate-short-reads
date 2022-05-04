#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { art_illumina as simulate_reads } from './modules/simulate_reads.nf'
include { simulate_contaminant_reads } from './modules/simulate_reads.nf'
include { downsample_simulated_reads } from './modules/simulate_reads.nf'
include { downsample_contaminant_reads } from './modules/simulate_reads.nf'
include { introduce_contaminants } from './modules/simulate_reads.nf'
include { fastp } from './modules/simulate_reads.nf'
include { bwa_align } from './modules/simulate_reads.nf'
include { downsample_regions } from './modules/simulate_reads.nf'
include { qualimap_bamqc } from './modules/simulate_reads.nf'
include { qualimap_bamqc_genome_results_to_csv } from './modules/simulate_reads.nf'
include { samtools_stats } from './modules/simulate_reads.nf'


workflow {
  ch_assemblies = Channel.fromPath( params.assembly_search_path ).map{ it -> [it.baseName, it] }.unique{ it -> it[0] }
  
  ch_replicates = Channel.fromList([1..params.replicates][0])

  if (params.depths_file != 'NO_FILE') {
    ch_depths = Channel.fromPath(params.depths_file).splitCsv(header: false)
  } else {
    ch_depths = Channel.of(params.depth)
  }

  if (params.downsample_regions_file != 'NO_FILE') {
    ch_downsample_regions = Channel.fromPath(params.downsample_regions_file).splitCsv(header: true).map{ it -> [it['ID'], it['DOWNSAMPLE_REGIONS']] }
  }

  if (params.contaminants != 'NO_FILE') {
    ch_contaminants = Channel.fromPath(params.contaminants).splitCsv(header: true).map{ it -> [it['ID'], it['ASSEMBLY'], it['PROPORTION']] }
    simulate_contaminant_reads(ch_contaminants)
    simulate_reads(ch_assemblies.combine(ch_depths).combine(ch_replicates))
    ch_proportion_uncontaminated = ch_contaminants.map{ it -> Float.parseFloat(it[2]) }.reduce{ x, y -> (x + y).round(6) }
    downsample_simulated_reads(simulate_reads.out.reads.combine(ch_proportion_uncontaminated))
    downsample_contaminant_reads(simulate_contaminant_reads.out.combine(simulate_reads.out.reads))
    introduce_contaminants(downsample_simulated_reads.out.join(downsample_contaminant_reads.out.uncompressed_reads.groupTuple(by: [0, 1]), by: [0, 1]))
    ch_reads = introduce_contaminants.out.reads
  } else {
    simulate_reads(ch_assemblies.combine(ch_depths).combine(ch_replicates))
    ch_reads = simulate_reads.out.reads
  }

  main:
    fastp(ch_reads)
    bwa_align(ch_assemblies.cross(ch_reads).map{ it -> [it[1][0], it[1][1], it[1][2], it[1][3], it[0][1]] })
    if (params.downsample_regions_file != 'NO_FILE') {
      downsample_regions(bwa_align.out.join(ch_downsample_regions))
      ch_aligned = downsample_regions.out
    } else {
      ch_aligned = bwa.align.out
    }
    samtools_stats(ch_aligned)
    qualimap_bamqc(ch_aligned)
    qualimap_bamqc_genome_results_to_csv(qualimap_bamqc.out.genome_results)
}
