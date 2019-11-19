#! /usr/bin/python3

import argparse
from Bio import SeqIO
import numpy as np
import re
import sys


def count_filtered(
    fastq,
    samples,
    kmers,
    filter_matrix,
    kmer_len,
    out_prefix,
    suffix
):
    '''
    Filters the FASTQ by filter matrix, writes to out FASTQ, and returns stats.
    '''
    name_convert = re.sub(".f[ast]*q", suffix, fastq)
    samples_index = samples.index(name_convert)
    #print(samples)
    #print(filter_matrix[0])
    forbidden = set([
        kmers[i]
        for i, x in enumerate(filter_matrix[:, samples_index])
        if x
    ])
    total = 0
    failed = 0
    out_records = []
    #out_name = out_prefix + fastq
    out_name = re.sub(".f[ast]*q", ".kmer_filtered.fq", fastq)
    for seq_record in SeqIO.parse(fastq, 'fastq'):
        seq = str(seq_record.seq)
        total += 1
        passes_filter = True
        for i in range(len(seq) - kmer_len + 1):
            kmer = seq[i: i + kmer_len]
            if kmer in forbidden:
                failed += 1
                passes_filter = False
                break
        if passes_filter:
            out_records.append(seq_record)
    SeqIO.write(out_records, out_name, 'fastq')
    return failed, total


def read_filter_tsv(filename):
    '''Parses TSV to return a matrix.'''
    filter_matrix = []
    kmer_vector = []
    with open(filename, 'r') as handle:
        samples = handle.readline().strip('\n').split('\t')
        for line in handle:
            cols = line.strip('\n').split('\t')
            kmer = cols[0]
            filter_vector = [bool(int(i)) for i in cols[1 : ]]
            kmer_vector.append(kmer)
            filter_matrix.append(filter_vector)
    return samples, kmer_vector, np.array(filter_matrix)
