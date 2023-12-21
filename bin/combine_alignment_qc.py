#!/usr/bin/env python3

import argparse
import csv
import json
import sys

def parse_samtools_stats_summary(samtools_stats_summary_file):
    samtools_stats_summary_data = []
    with open(samtools_stats_summary_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            record = {}
            for k, v in row.items():
                new_k = k + '__samtools'
                record[new_k] = v
            samtools_stats_summary_data.append(record)

    first_record = samtools_stats_summary_data[0]
    
    return first_record


def parse_qualimap_bamqc_genome_results(qualimap_bamqc_genome_results_file):
    qualimap_bamqc_genome_results_data = []
    with open(qualimap_bamqc_genome_results_file) as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            record = {}
            for k, v in row.items():
                new_k = k + '__qualimap'
                record[new_k] = v
            qualimap_bamqc_genome_results_data.append(record)

    first_record = qualimap_bamqc_genome_results_data[0]
    
    return first_record


def main(args):
    samtools_stats_summary_data = parse_samtools_stats_summary(args.samtools_stats_summary)
    qualimap_bamqc_genome_results_data = parse_qualimap_bamqc_genome_results(args.qualimap_bamqc_genome_results)

    samtools_fields = [
        'reads_mapped',
        'reads_mapped_and_paired',
        'reads_unmapped',
        'reads_properly_paired',
        'reads_paired',
        'reads_duplicated',
        'error_rate',
        'non-primary_alignments',
        'supplementary_alignments',
        'average_length',
        'average_first_fragment_length',
        'average_last_fragment_length',
        'insert_size_average',
        'insert_size_standard_deviation',
    ]

    qualimap_fields = [
        'mean_depth_coverage',
        'stdev_depth_coverage',
    ]

    selected_samtools_data = { k: samtools_stats_summary_data[k + '__samtools'] for k in samtools_fields }
    selected_qualimap_data = { k: qualimap_bamqc_genome_results_data[ k + '__qualimap'] for k in qualimap_fields }
    combined_data = selected_samtools_data.copy()
    combined_data.update(selected_qualimap_data)

    output_fields = []
    
    if args.sample_id:
        combined_data['sample_id'] = args.sample_id
        output_fields.append('sample_id')
    if args.read_type:
        combined_data['read_type'] = args.read_type
        output_fields.append('read_type')

    output_fields += samtools_fields
    output_fields += qualimap_fields
  
    writer = csv.DictWriter(sys.stdout, dialect='unix', fieldnames=output_fields, quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    writer.writerow(combined_data)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--samtools-stats-summary', required=True)
    parser.add_argument('--qualimap-bamqc-genome-results', required=True)
    parser.add_argument('--sample-id')
    parser.add_argument('--read-type')
    args = parser.parse_args()
    main(args)
