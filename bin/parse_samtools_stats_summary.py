#!/usr/bin/env python3

import argparse
import csv
import json
import sys


def parse_samtools_stats_summary(samtools_stats_summary_file):
    output_data = {}
    int_fields = [
        'raw_total_sequences',
        'filtered_sequences',
        'sequences',
        'first_fragments',
        'last_fragments',
        'reads_mapped',
        'reads_mapped_and_paired',
        'reads_unmapped',
        'reads_properly_paired',
        'reads_paired',
        'reads_duplicated',
        'reads_MQ0',
        'reads_QC_failed',
        'non-primary_alignments',
        'supplementary_alignments',
        'total_length',
        'total_first_fragment_length',
        'total_last_fragment_length',
        'bases_mapped',
        'bases_mapped_cigar',
        'bases_trimmed',
        'bases_duplicated',
        'mismatches',
        'maximum_length',
        'maximum_first_fragment_length',
        'maximum_last_fragment_length',
        'inward_oriented_pairs',
        'outward_oriented_pairs',
        'pairs_with_other_orientation',
        'pairs_on_different_chromosomes',
    ]
    float_fields = [
        'error_rate',
        'average_length',
        'average_first_fragment_length',
        'average_last_fragment_length',
        'average_quality',
        'insert_size_average',
        'insert_size_standard_deviation',
    ]
    with open(samtools_stats_summary_file, 'r') as f:
        for line in f:
            line = line.strip()
            line_split = line.split('\t')
            key = line_split[0].strip().rstrip(':').replace(' ', '_').replace('(', '').replace(')', '')
            if key == '1st_fragments':
                key = 'first_fragments'
            if key.endswith('%'):
                key = key.replace('_%', '')
            elif key in int_fields:
                try:
                    output_data[key] = int(line_split[1])
                except ValueError:
                    output_data[key] = None
            elif key == 'is_sorted':
                output_data[key] = line_split[1] == '1'
            elif key in float_fields:
                try:
                    output_data[key] = float(line_split[1])
                except ValueError:
                    output_data[key] = None
            else:
                output_data[key] = line_split[1]

    return output_data


def main(args):
    samtools_stats_summary_data = parse_samtools_stats_summary(args.input)

    if args.sample_id:
        samtools_stats_summary_data['sample_id'] = args.sample_id

    output_fields = [
        'raw_total_sequences',
        'filtered_sequences',
        'first_fragments',
        'last_fragments',
        'reads_mapped',
        'reads_mapped_and_paired',
        'reads_unmapped',
        'reads_properly_paired',
        'reads_paired',
        'reads_duplicated',
        'reads_MQ0',
        'reads_QC_failed',
        'non-primary_alignments',
        'supplementary_alignments',
        'total_length',
        'total_first_fragment_length',
        'total_last_fragment_length',
        'bases_mapped',
        'bases_mapped_cigar',
        'bases_trimmed',
        'bases_duplicated',
        'mismatches',
        'error_rate',
        'average_length',
        'average_first_fragment_length',
        'average_last_fragment_length',
        'average_quality',
        'insert_size_average',
        'insert_size_standard_deviation',
        'inward_oriented_pairs',
        'outward_oriented_pairs',
        'pairs_with_other_orientation',
        'pairs_on_different_chromosomes',
    ]
    if args.sample_id:
        output_fields = ['sample_id'] + output_fields

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    writer.writerow(samtools_stats_summary_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse samtools stats summary file')
    parser.add_argument('-i', '--input', type=str, required=True, help='samtools stats summary file')
    parser.add_argument('-s', '--sample-id', type=str, required=True, help='sample id')
    args = parser.parse_args()
    main(args)
