#!/usr/bin/env python

import argparse
import csv
import json
import sys

def main(args):
    with open(args.qualimap_bamqc_genome_results, 'r') as f:
        output_data = {
            'sample_id': args.sample_id,
            'read_type': args.read_type,
        }
        for line in f:
            line = line.strip()
            if line.startswith('number of reads'):
                number_of_reads = line.split('=')[1].strip().replace(',', '')
                output_data['num_reads'] = int(number_of_reads)
            if line.startswith('number of mapped reads'):
                num_mapped_reads = line.split('=')[1].strip().split(' ')[0].replace(',', '')
                output_data['num_mapped_reads'] = int(num_mapped_reads)
                percent_mapped_reads = line.split('=')[1].strip().split(' ')[1].strip().replace('(', '').replace(')', '').replace('%', '')
                output_data['percent_mapped_reads'] = round(float(percent_mapped_reads), 2)
            if line.startswith('number of secondary alignments'):
                num_secondary_alignments = int(line.split('=')[1].strip().replace(',', ''))
                output_data['num_secondary_alignments'] = num_secondary_alignments
            if line.startswith('duplication rate'):
                duplication_rate = line.split('=')[1].strip().replace('%', '')
                output_data['duplication_rate_percent'] = round(float(duplication_rate), 2)
            if line.startswith('mean coverageData'):
                mean_coverage = line.split('=')[1].strip().strip('X').replace(',', '')
                output_data['mean_depth_coverage'] = round(float(mean_coverage), 2)
            if line.startswith('std coverageData'):
                stdev_coverage = line.split('=')[1].strip().strip('X').replace(',', '')
                output_data['stdev_depth_coverage'] = round(float(stdev_coverage), 2)
            if line.startswith('mean mapping quality'):
                mean_mapping_quality = line.split('=')[1].strip()
                output_data['mean_mapping_quality'] = round(float(mean_mapping_quality), 2)
            if line.startswith('general error rate'):
                general_error_rate = line.split('=')[1].strip()
                output_data['error_rate'] = round(float(general_error_rate), 2)
            if line.startswith('number of mismatches'):
                number_of_mismatches = line.split('=')[1].strip().replace(',', '')
                output_data['number_of_mismatches'] = int(number_of_mismatches)
            if line.startswith('number of insertions'):
                number_of_insertions = line.split('=')[1].strip().replace(',', '')
                output_data['number_of_insertions'] = int(number_of_insertions)
            if line.startswith('mapped reads with insertion percentage'):
                mapped_reads_with_insertion_percentage = line.split('=')[1].strip().replace('%', '')
                output_data['mapped_reads_with_insertion_percentage'] = round(float(mapped_reads_with_insertion_percentage), 2)
            if line.startswith('number of deletions'):
                number_of_deletions = line.split('=')[1].strip().replace(',', '')
                output_data['number_of_deletions'] = int(number_of_deletions)
            if line.startswith('mapped reads with deletion percentage'):
                mapped_reads_with_deletion_percentage = line.split('=')[1].strip().replace('%', '')
                output_data['mapped_reads_with_deletion_percentage'] = round(float(mapped_reads_with_deletion_percentage), 2)
            if 'reference with a coverageData >= 5X' in line:
                proportion_genome_covered_over_5x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_5x'] = round(proportion_genome_covered_over_5x, 4)
            if 'reference with a coverageData >= 10X' in line:
                proportion_genome_covered_over_10x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_10x'] = round(proportion_genome_covered_over_10x, 4)
            if 'reference with a coverageData >= 20X' in line:
                proportion_genome_covered_over_20x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_20x'] = round(proportion_genome_covered_over_20x, 4)
            if 'reference with a coverageData >= 30X' in line:
                proportion_genome_covered_over_30x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_30x'] = round(proportion_genome_covered_over_30x, 4)
            if 'reference with a coverageData >= 40X' in line:
                proportion_genome_covered_over_40x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_40x'] = round(proportion_genome_covered_over_40x, 4)
            if 'reference with a coverageData >= 50X' in line:
                proportion_genome_covered_over_50x = float(line.split(' ')[3].strip('%')) / 100
                output_data['proportion_genome_covered_over_50x'] = round(proportion_genome_covered_over_50x, 4)
            

    output_fieldnames = [
        'sample_id',
        'read_type',
        'mean_depth_coverage',
        'stdev_depth_coverage',
        'num_reads',
        'num_mapped_reads',
        'percent_mapped_reads',
        'mean_mapping_quality',
        'error_rate',
        'number_of_mismatches',
        'number_of_insertions',
        'mapped_reads_with_insertion_percentage',
        'number_of_deletions',
        'mapped_reads_with_deletion_percentage',
        'num_secondary_alignments',
        'duplication_rate_percent',
        'proportion_genome_covered_over_5x',
        'proportion_genome_covered_over_10x',
        'proportion_genome_covered_over_20x',
        'proportion_genome_covered_over_30x',
        'proportion_genome_covered_over_40x',
        'proportion_genome_covered_over_50x',
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', extrasaction='ignore', quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerow(output_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('qualimap_bamqc_genome_results')
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-t', '--read-type', choices=['short', 'long'], default='short')
    args = parser.parse_args()
    main(args)
