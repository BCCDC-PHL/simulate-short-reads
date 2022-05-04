#!/usr/bin/env python

import argparse
import csv
import json
import random
import re
import sys

import pysam

from collections import defaultdict
from functools import reduce


def get_idxstats(bam_path):
    idxstats = {
        "ref_name": "",
        "seq_length": 0,
        "num_mapped_segments": 0,
        "num_unmapped_segments": 0,
    }

    idxstats_out = pysam.idxstats(bam_path).strip().split('\t')
    idxstats_out[3] = idxstats_out[3].split('\n')[0]

    idxstats['ref_name'] = idxstats_out[0]
    idxstats['seq_length'] = int(idxstats_out[1])
    idxstats['num_mapped_segments'] = int(idxstats_out[2])
    idxstats['num_unmapped_segments'] = int(idxstats_out[3])

    return idxstats



def read_pair_generator(bam, region_string=None):
    """
    from: https://www.biostars.org/p/306041/#332022
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_downsampling_bed(downsampling_bed_path):
    fieldnames = [
        'chrom',
        'chromStart',
        'chromEnd',
        'name',
        'score',
    ]
    downsampling_bed_by_chrom = {}
    with open(downsampling_bed_path, 'r') as f:
        reader = csv.DictReader(f, fieldnames=fieldnames, dialect='excel-tab')
        for row in reader:
            row['chromStart'] = int(row['chromStart'])
            row['chromEnd'] = int(row['chromEnd'])
            row['score'] = float(row['score'])
            if row['chrom'] not in downsampling_bed_by_chrom:
                downsampling_bed_by_chrom[row['chrom']] = [row]
            else:
                downsampling_bed_by_chrom[row['chrom']].append(row)

    return downsampling_bed_by_chrom


def write_segments(segment, mate_segment, outfile):
    outfile.write(segment)
    outfile.write(mate_segment)

    return True


def maybe_write_segments(segment, mate_segment, outfile, probability_of_writing):
    random_num = random.uniform(0, 1)
    will_write = random_num < probability_of_writing
    if will_write:
        outfile.write(segment)
        outfile.write(mate_segment)

    return will_write


def point_in_region(point, region_min, region_max):
    in_region = False
    if (point > region_min) and (point < region_max):
        in_region = True
    return in_region


def main(args):
    """
    """

    idxstats = get_idxstats(args.bam)
    total_reads = idxstats['num_mapped_segments'] + idxstats['num_unmapped_segments']

    # open the primer scheme and get the pools
    downsampling_bed = read_downsampling_bed(args.bed)

    infile = pysam.AlignmentFile(args.bam, "rb")
    bam_header = infile.header.copy().to_dict()

    outfile = pysam.AlignmentFile("-", "wh", header=bam_header)

    for segment, mate_segment in read_pair_generator(infile):
        segments_written = False
        segments_in_a_region = False
        for chrom in downsampling_bed:
            if (segment.reference_name == chrom) and (mate_segment.reference_name == chrom):
                for region in downsampling_bed[chrom]:
                    segment_in_region = point_in_region(segment.reference_start, region['chromStart'], region['chromEnd']) or point_in_region(segment.reference_end, region['chromStart'], region['chromEnd'])
                    if segment_in_region:
                        segments_in_a_region = True
                        segments_written = maybe_write_segments(segment, mate_segment, outfile, region['score'])

                if not segments_in_a_region:
                    write_segments(segment, mate_segment, outfile)

    # close up the file handles
    infile.close()
    outfile.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Downsample alignments from an amplicon scheme.')
    parser.add_argument('bam', help='bam file containing the alignment')
    parser.add_argument('--bed', help='BED file containing the downsampling regions')
    parser.add_argument('--verbose', action='store_true', help='Debug mode')
    args = parser.parse_args()
    main(args)
