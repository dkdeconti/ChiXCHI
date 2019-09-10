#! /usr/bin/python3

import argparse
import sys
import ChiXCHI.common.kmer_filter as kmer_filter
import ChiXCHI.common.filter_reads as flt_reads


def identify_kmers(args):
    '''
    Main function for handling args for identify subparser.
    '''
    kmers, kmer_matrix = kmer_filter.normalize_kmers(
        kmer_filter.create_kmer_matrix(args.kmers)
    )
    fit_vector = kmer_filter.test_for_fitness(kmer_matrix, args.degrees)
    filter_vector = kmer_filter.cluster_kmers(kmer_matrix, fit_vector, args.threads)
    kmer_filter.write_matrix(filter_vector, kmers, args.kmers)


def filter_contaminants(args):
    '''
    Main function for handling args for the filtering of FASTQs given the
    previously identified contaminant k-mers.
    '''
    samples, kmers, filter_matrix = flt_reads.read_filter_tsv(args.filter_tsv)
    results = []
    for fastq in args.fastqs:
        results.append(flt_reads.count_filtered(
            fastq,
            samples,
            kmers,
            filter_matrix,
            args.length,
            args.output,
            args.suffix
        ))
    results = [failed / float(total) for failed, total in results]


def main():
    '''
    Handles argument parsing and dispatch of methods.
    '''
    parser = argparse.ArgumentParser(prog='ChiXCHI')
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')
    # Create the parser for identifying k-mers clusters
    parser_identify = subparsers.add_parser(
        'identify',
        help='identify help'
    )
    parser_identify.add_argument(
        '-t',
        '--threads',
        metavar='THREADS',
        type=int,
        default=1,
        help='number of threads to use'
    )
    parser_identify.add_argument(
        '-d',
        '--degrees',
        metavar='DEGREES_OF_FREEDOM',
        type=int,
        default=0,
        help='degrees of freedom for chi-squared test'
    )
    parser_identify.add_argument(
        'kmers',
        metavar='KMER_TSV',
        type=str,
        nargs='+',
        help="k-mer count TSV"
    )
    # Create the parser for filtering reads
    # ToDo: Fix arguments for filter subparser
    parser_filter = subparsers.add_parser(
        'filter',
        help='filter help'
    )
    parser_filter.add_argument(
        '-s',
        '--suffix',
        metavar='SUFFIX',
        default='.kmers.tsv',
        help='suffix for matrix header; default: .kmers.tsv'
    )
    parser_filter.add_argument(
        '-l',
        '--length',
        metavar='KMER_LENGTH',
        type=int, default=31,
        help='length of k-mer; default matches BFCounter default: 31'
    )
    parser_filter.add_argument(
        '-o',
        '--output',
        metavar='OUPUT_PREFIX',
        default='filtered_',
        help='output prefix for FASTA filename'
    )
    parser_filter.add_argument(
        'tsv',
        metavar='TSV',
        help="filter TSV matrix for keywords from sub-command \'identify\'"
    )
    parser_filter.add_argument(
        'fastqs',
        metavar='FASTQ',
        nargs='+',
        help='FASTQ file(s) to be filtered'
    )
    #
    args = parser.parse_args()
    dispatch = {
        'identify': identify_kmers,
        'filter': filter_contaminants
    }
    dispatch[args.command](args)
    


if __name__ == "__main__":
    main()
