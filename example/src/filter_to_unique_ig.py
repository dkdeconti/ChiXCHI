#! /usr/bin/env python
"""Identifies the unique sequences in ig_simulator's final repertoire FASTA."""


import argparse
import sys
from Bio import SeqIO


def uniqify(filename):
    """Parses FASTA and uniqifies antibodies and returns seq records list."""
    seq_dict = {}
    try:
        for seq_record in SeqIO.parse(filename, "fasta"):
            new_id = "_".join(seq_record.id.split('_')[:2])
            seq_record.id = new_id
            seq_record.description = ""
            if seq_record.id not in seq_dict:
                seq_dict[seq_record.id] = seq_record
    except:
        sys.stderr.write("There was an issue parsing the FASTA file.")
        sys.exit()
    return seq_dict.values()


def write_seq_records(seq_records, output_filename):
    """Writes seq records as FASTA to stdout or output filename."""
    if output_filename:
        SeqIO.write(seq_records, output_filename, "fasta")
    else:
        SeqIO.write(seq_records, sys.stdout, "fasta")


def main():
    """Arg parsing and central dispatch"""
    # arg parsing
    parser = argparse.ArgumentParser(description="Uniquify IG repertoire.")
    parser.add_argument("fasta", metavar="FASTA",
                        help="final repertoire FASTA from ig_simulator")
    parser.add_argument("-o", "--output", metavar="FILENAME",
                        help="output FASTA filename")
    args = parser.parse_args()
    # central dispatch
    seq_records = uniqify(args.fasta)
    write_seq_records(seq_records, args.output)


if __name__ == "__main__":
    main()
