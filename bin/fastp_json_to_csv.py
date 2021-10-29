#!/usr/bin/env python

import argparse
import json

def main(args):
    with open(args.fastp_json, 'r') as f:
        fastp_report = json.load(f)

    total_reads_before_filtering = fastp_report['summary']['before_filtering']['total_reads']
    total_reads_after_filtering = fastp_report['summary']['after_filtering']['total_reads']
    total_bases_before_filtering = fastp_report['summary']['before_filtering']['total_bases']
    total_bases_after_filtering = fastp_report['summary']['after_filtering']['total_bases']
    q20_bases_before_filtering = fastp_report['summary']['before_filtering']['q20_bases']
    q20_bases_after_filtering = fastp_report['summary']['after_filtering']['q20_bases']
    q30_bases_before_filtering = fastp_report['summary']['before_filtering']['q30_bases']
    q30_bases_after_filtering = fastp_report['summary']['after_filtering']['q30_bases']
    adapter_trimmed_reads = fastp_report['adapter_cutting']['adapter_trimmed_reads']
    adapter_trimmed_bases = fastp_report['adapter_cutting']['adapter_trimmed_bases']

    output_fields = [
        'total_reads_before_filtering',
        'total_reads_after_filtering',
        'total_bases_before_filtering',
        'total_bases_after_filtering',
        'q20_bases_before_filtering',
        'q20_bases_after_filtering',
        'q30_bases_before_filtering',
        'q30_bases_after_filtering',
        'adapter_trimmed_reads',
        'adapter_trimmed_bases',
    ]

    output_data = []
    if args.sample_id:
        output_fields = ['sample_id'] + output_fields
        output_data = [args.sample_id]

    print(",".join(output_fields))
    output_data = output_data + [
        total_reads_before_filtering,
        total_reads_after_filtering,
        total_bases_before_filtering,
        total_bases_after_filtering,
        q20_bases_before_filtering,
        q20_bases_after_filtering,
        q30_bases_before_filtering,
        q30_bases_after_filtering,
        adapter_trimmed_reads,
        adapter_trimmed_bases,
    ]
    print(",".join(map(str, output_data)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastp_json')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
