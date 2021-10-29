#!/usr/bin/env python

import argparse
import csv
import json
import sys

def main(args):
    with open(args.qualimap_bamqc_genome_results, 'r') as f:
        output_data = { 'sample_id': args.sample_id }
        for line in f:
            line = line.strip()
            if line.startswith('median insert size'):
                median_insert_size = line.split('=')[1].strip().replace(',', '')
                output_data['median_insert_size'] = int(median_insert_size)
            if line.startswith('mean coverageData'):
                mean_coverage = line.split('=')[1].strip().strip('X')
                output_data['mean_coverage'] = round(float(mean_coverage), 2)
            if line.startswith('std coverageData'):
                stdev_coverage = line.split('=')[1].strip().strip('X')
                output_data['stdev_coverage'] = round(float(stdev_coverage), 2)
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
        'median_insert_size',
        'mean_coverage',
        'stdev_coverage',
        'proportion_genome_covered_over_5x',
        'proportion_genome_covered_over_10x',
        'proportion_genome_covered_over_20x',
        'proportion_genome_covered_over_30x',
        'proportion_genome_covered_over_40x',
        'proportion_genome_covered_over_50x',
    ]

    csv.register_dialect('unix-csv', delimiter=',', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix-csv')
    writer.writeheader()
    writer.writerow(output_data)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('qualimap_bamqc_genome_results')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
