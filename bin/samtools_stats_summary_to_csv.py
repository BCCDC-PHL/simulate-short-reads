#!/usr/bin/env python

import argparse
import csv
import json
import sys

def parse_samtools_stats_summary(samtools_stats_summary_path):
    s = { 'sample_id': args.sample_id }
    with open(samtools_stats_summary_path, 'r') as f:   
        for line in f:
            line = line.strip()
            key = line.split('\t')[0].lower().strip(':').replace(' ', '_').replace('(','').replace(')','').replace('%', 'percent')
            value = line.split('\t')[1]
            s[key] = value

    return s

def main(args):
    samtools_stats_summary = parse_samtools_stats_summary(args.samtools_stats_summary)
    

    output_fieldnames = [
        'sample_id',
        'raw_total_sequences',
        'bases_trimmed',
        'error_rate',
        'insert_size_standard_deviation',
        'non-primary_alignments',
        'sequences',
        'maximum_last_fragment_length',
        'reads_mapped_and_paired',
        'average_quality',
        'insert_size_average',
        'reads_mapped',
        'reads_unmapped',
        'reads_qc_failed',
        'reads_paired',
        'average_last_fragment_length',
        '1st_fragments',
        'reads_mq0',
        'average_length',
        'supplementary_alignments',
        'bases_mapped',
        'total_length',
        'bases_mapped_cigar',
        'bases_duplicated',
        'average_first_fragment_length',
        'pairs_with_other_orientation',
        'mismatches',
        'total_last_fragment_length',
        'reads_properly_paired',
        'inward_oriented_pairs',
        'filtered_sequences',
        'outward_oriented_pairs',
        'percentage_of_properly_paired_reads_percent',
        'total_first_fragment_length',
        'pairs_on_different_chromosomes',
        'maximum_length',
        'last_fragments',
        'is_sorted',
        'maximum_first_fragment_length',
        'reads_duplicated',
    ]

    csv.register_dialect('unix-csv', delimiter=',', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix-csv')
    writer.writeheader()
    writer.writerow(samtools_stats_summary)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('samtools_stats_summary')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
