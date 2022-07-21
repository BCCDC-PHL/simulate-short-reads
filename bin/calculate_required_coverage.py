#!/usr/bin/env python3

import argparse


def main(args):
    num_contaminant_reads_required = int(args.num_simulated_reads * args.contaminant_proportion)
    
    depth_coverage_required = int(num_contaminant_reads_required * args.read_length / args.contaminant_genome_size)

    depth_coverage_to_generate = int(depth_coverage_required * args.excess_reads_proportion)

    if depth_coverage_to_generate < 1:
        depth_coverage_to_generate = 1

    print(depth_coverage_to_generate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--num-simulated-reads', type=int)
    parser.add_argument('--read-length', type=int)
    parser.add_argument('--contaminant-genome-size', type=int)
    parser.add_argument('--contaminant-proportion', type=float)
    parser.add_argument('--excess-reads-proportion', type=float, default=1.5)
    args = parser.parse_args()
    main(args)
