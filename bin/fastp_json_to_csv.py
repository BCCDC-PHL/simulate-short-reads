#!/usr/bin/env python

import argparse
import json

def main(args):
    with open(args.fastp_json, 'r') as f:
        fastp_report = json.load(f)

    total_reads = fastp_report['summary']['before_filtering']['total_reads']
    total_bases = fastp_report['summary']['before_filtering']['total_bases']
    q20_bases = fastp_report['summary']['before_filtering']['q20_bases']
    q30_bases = fastp_report['summary']['before_filtering']['q30_bases']
    q20_rate = fastp_report['summary']['before_filtering']['q20_rate']
    q30_rate = fastp_report['summary']['before_filtering']['q30_rate']
        
    output_fields = [
        'total_reads',
        'total_bases',
        'q20_bases',
        'q30_bases',
        'q20_rate',
        'q30_rate',
    ]

    output_data = []
    if args.sample_id:
        output_fields = ['sample_id'] + output_fields
        output_data = [args.sample_id]

    print(",".join(output_fields))
    output_data = output_data + [
        total_reads,
        total_bases,
        q20_bases,
        q30_bases,
        q20_rate,
        q30_rate,
    ]
    print(",".join(map(str, output_data)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastp_json')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
