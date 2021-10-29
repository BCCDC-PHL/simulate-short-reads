#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { simulate_reads } from './modules/simulate_reads.nf'
include { fastp } from './modules/simulate_reads.nf'
include { fastp_json_to_csv } from './modules/simulate_reads.nf'
include { bwa_align } from './modules/simulate_reads.nf'
include { qualimap_bamqc } from './modules/simulate_reads.nf'
include { qualimap_bamqc_genome_results_to_csv } from './modules/simulate_reads.nf'


workflow {
  ch_assemblies = Channel.fromPath( params.assembly_search_path ).map{ it -> [it.baseName, it] }.unique{ it -> it[0] }
  ch_depths = Channel.fromList([25, 50, 75, 100])
  ch_replicates = Channel.fromList([1..params.replicates][0])

  main:
    simulate_reads(ch_assemblies.combine(ch_depths).combine(ch_replicates))
    fastp(simulate_reads.out.reads)
    fastp_json_to_csv(fastp.out.json)
    bwa_align(ch_assemblies.cross(simulate_reads.out.reads).map{ it -> [it[1][0], it[1][1], it[1][2], it[1][3], it[0][1]] })
    qualimap_bamqc(bwa_align.out)
    qualimap_bamqc_genome_results_to_csv(qualimap_bamqc.out.genome_results)
}
