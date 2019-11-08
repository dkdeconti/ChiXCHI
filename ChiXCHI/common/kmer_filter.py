#! /usr/bin/python3

import argparse
import collections
import numpy as np
import sklearn.cluster
import sklearn.preprocessing
import scipy
import sys
from multiprocessing import Pool
from statsmodels.sandbox.stats.multicomp import multipletests


def cluster_kmers(kmer_freqs, padj, threads):
    '''
    Cluster kmers not filtered by chi-squared goodness of fit to return
    a filter matrix.
    '''
    filter_matrix = [[0]*len(v) for v in kmer_freqs]
    pool = Pool(threads)
    to_cluster = [(i, v) for i, v in enumerate(kmer_freqs) if padj[i]]
    # ToDo: figure out where missing code went...
    clusters = pool.map(cluster_vector, to_cluster)
    pool.close()
    for i, v in clusters:
        filter_matrix[i] = v
    return filter_matrix


def cluster_vector(pooled_args):
    '''
    Performs clustering on k-mer frequences and returns vector with labeled
    failed cluster.
    '''
    i, v = pooled_args
    ap = sklearn.cluster.AffinityPropagation().fit(v.reshape(-1, 1))
    fail_cluster = ap.labels_[np.argmin(v)]
    return i, [clust == fail_cluster for clust in ap.labels_]


def create_kmer_matrix(filenames):
    '''Parses TSV k-mer count file to return a count matrix as dict.'''
    kmer_matrix = collections.defaultdict(lambda: [0]*len(filenames))
    for idx, filename in enumerate(filenames):
        with open(filename, 'r') as handle:
            for line in handle:
                cols = line.strip('\n').split()
                kmer, count = cols
                try:
                    count = int(count)
                except ValueError:
                    sys.stderr.write('ValueError typecasting k-mer count.\n')
                    sys.exit()
                kmer_matrix[kmer][idx] = count
    return kmer_matrix


def normalize_kmers(nonnormalized_matrix, max_norm_values=100000):
    '''Normalizes matrix to a minmax scale'''
    kmers = nonnormalized_matrix.keys()
    counts = np.array([nonnormalized_matrix[k] for k in kmers])
    counts_minmax = sklearn.preprocessing.minmax_scale(
        counts,
        feature_range=(0, max_norm_values)
    )
    return list(kmers), counts_minmax


def test_for_fitness(count_matrix, df):
    '''
    Parses vecotr for 1st pass with chi-squared on normalized count frequences
    to return a pre-filtered matrix of non-contaminant k-mers.
    '''
    pvalues = scipy.stats.chisquare(count_matrix, axis=1, ddof=df)[1]
    return multipletests(pvalues, method='fdr_bh')[0]


def write_matrix(filter_vector, kmers, tsv_filenames):
    '''Writes matrix to stdout and prepends kmer.'''
    sys.stdout.write('\t'.join(tsv_filenames) + '\n')
    for i, v in enumerate(filter_vector):
        try:
            output_v = [str(int(b)) for b in v]
        except ValueError:
            e = "ValueError typecasting filter vector in write_matrix.\n"
            sys.stderr.write(e)
            sys.exit()
        sys.stdout.write('\t'.join([kmers[i]] + output_v) + '\n')