# ChiXCHI
  Software for identifying cross-contaminant samples in high throughput sequencing.

## Overview

ChiXCHI was designed primarily for identifying cross-sample and chimeric read contaminants for PCR based amplicon next gen sequencing.

## Requirements

* Python3
* 

### Optional software

* bbmerge <sup>[1, ](#bbmerge)</sup><sup>[2](#bbmergsite)</sup>
* BFCounter <sup>[3, ](#bfcounter)</sup><sup>[4](#bfcountergithub)</sup>

## Pre-processing

ChiXCHI requires some pre-processing steps to for appropriate input data:

1. FASTQ
2. Two column text file of k-mer counts 

It is also recommended, but not required, that reads be merged prior to analysis. There is currently no built-in methods for handling filtering of read pairs. bbmerge <sup>[1, ](#bbmerge)</sup><sup>[2](#bbmergsite)</sup> is one recommended tool for merging paired FASTQs. BFCounter <sup>[3, ](#bfcounter)</sup><sup>[4](#bfcountergithub)</sup> is a recommended tool for producing the necessary table for k-mer counts.

## Usage

```bash

```

## Example usage with synthetic data

Here I outline an example of running the software with a synthetically created IG dataset. This setup will show the ability to identify inter-well contamination rates. Running this example will require:

* python2.7
* bbmerge <sup>[1, ](#bbmerge)</sup><sup>[2](#bbmergsite)</sup>
* BFCounter <sup>[3, ](#bfcounter)</sup><sup>[4](#bfcountergithub)</sup>
* IgSimulator <sup>[5, ](#igsimulator)</sup><sup>[6](#igsimulatorgithub)</sup>
* ART <sup>[7, ](#art)</sup><sup>[8](#artdl)</sup>
* seqtk <sup>[9](#seqtk)</sup>

All other scripts can be found in example/src.

### Create the repertoire of simulated IGs

```bash
# Create a simulated repertoire of IGs.
python ig_simulator.py \
    --chain-type HC \
    --num-bases 1000 \
    --num-mutated 2000 \
    --repertoire-size 15000 \
    -o simulated_ig;

# Move the output repertoire into the current working directory
mv simulated_ig/final_repertoire.fasta .

# Filter the FASTA of IGs to only unique ones.
python example/src/filter_to_unique_ig.py \
    final_repertoire.fasta >> final_reportoire.unique.fasta
# Parameters for subsampling reads for a set of four randomly selected
# synthetic IGs.
num_wells=100;
min_subvars=90;
max_subvars=110;
num_subvars="range(${min_subvars}, ${max_subvars})";
cov_choices="[10000, 30000, 50000, 80000]";

# For the designated number of wells, we will assign that well a random
# number of unique IGs.
# 
# Then for each variant of IGs per well, we create a random coverage of
# the IG to synthetically sequence with ART. We merge the paired FASTQs
# with bbmerge and combine all the reads for that well into a single
# FASTQ.
for ((i=1; i <= $num_wells; i++)); do
    # Randomly select the number of sub-variants per well
    echo "#####################################################################"
    echo "Starting work on well ${i}..."
    echo "#####################################################################"
    num_vars=$(python -c "import random; print random.choice(${num_subvars});");
    ab_tracker="well_${i}\t";
    for ((j=$antibody_counter; j < antibody_counter + num_vars; j++)); do
        # Get single antibody
        tmp_antibodies="${ab_tracker}antibody_${j},";
        ab_tracker=$tmp_antibodies;
        echo "antibody_${j}" \
        | seqtk subseq final_repertoire.unique.fasta - \
        > tmp.ref.fa;
        # Randomly pick seq depth
        cov=$(python -c "import random; print random.choice(${cov_choices})");
        # Simulate reads with ART
        art_illumina \
            -ss MSv3 \
            -amp \
            -mp \
            -na \
            -l 250 \
            -f ${cov} \
            -i tmp.ref.fa \
            -o tmp.amp;
        # Merge paired mates with bbmerge
        bbmerge-auto.sh \
            in=tmp.amp1.fq \
            in2=tmp.amp2.fq \
            out=tmp.merged.fq \
            rem k=62 extend2=50 ecct;
        # Combine individual sub-variants into single well FASTQ
        cat tmp.merged.fq >> well_${i}.merged.fq;
    done;
    antibody_counter=$(($antibody_counter + $num_vars))
    echo -e "$ab_tracker" >> antibodies_tracker.tsv;
done;

# Parameters for setting the contamination rate between wells
baseline_depth=100000;
mu=$baseline_depth;
sigma=20000;
contam_rates="[x * 0.01 for x in range(1, 11)]";

# Now we contaminate the reads between wells, while also normalizing the
# library size for all wells. Library size is normalized around mu and sigma.
for well_fastq in $(ls well_*.merged.fq); do
    new_fastq=${well_fastq%.fq}.contam.fq;
    # Calculate the total depth.
    total_depth=$(python -c "import numpy; print numpy.random.normal(${mu}, ${sigma});");
    # Randomly select some contamination rate
    contam_rate=$(python -c "import random; print random.choice(${contam_rates});");
    # Identify the number of reads for contamination
    contam_depth=$(python -c "print int(${total_depth} * ${contam_rate});");
    # Identify the number of reads that are non-contamination
    real_depth=$(python -c "print int(${total_depth} * (1 - ${contam_rate}));");
    echo "$real_depth $contam_depth"
    # Subsample with seqtk the non-contamination reads
    seqtk sample $well_fastq $real_depth >> $new_fastq;
    # Subsample with seqtk the contamination reads from their respective
    # wells FASTQ. This is assumed random from the cohort of potential
    # contamination targets.
    cat $(ls well_*.merged.fq | grep -v ${well_fastq}) \
    | seqtk sample - $contam_depth >> $new_fastq;
done;
```

### Run ChiXCHI on the reads for 

```bash
# Count the k-mers with BFCounter
for fastq in $(ls well*.contam.fq); do
    ~/bin/BFCounter-master/BFCounter \
        count \
        -k 31 \
        -n 500000 \
        -o ${fastq%.fq}.kmers.bin \
        $fastq;
    ~/bin/BFCounter-master/BFCounter \
        dump \
        -k 31 \
        -i ${fastq%.fq}.kmers.bin \
        -o ${fastq%.fq}.kmers.tsv;
done;

# Identify the contaminant kmers
python ChiXCHI/identify_contaminants.py $(ls *.kmers.tsv) > filter_matrix.tsv;

# Filter the FASTQs for contaminant kmers
python ../filter_contaminants.py \
    -f filter_matrix.tsv \
    -l 31 \
    -o filtered_ \
    $(ls *.contam.fq);

# Calculate stats
echo -e "\taccuracy\tspecificity\tsensitivity\tprecision\tfdr" >> ../100_sample_4per.tsv
for fastq in well*.contam.fq; do
    python example/src/calculate_stats.py \
        -t antibodies_tracker.tsv \
        $fastq \
        filtered_${fastq} >> 100_sample_4per.tsv;
done;
```

## References

<a name="bbmerge">[1]</a> Bushnell B, Rood J, Singer E (2017) BBMerge – Accurate paired shotgun read merging via overlap. PLOS ONE 12(10): e0185056. https://doi.org/10.1371/journal.pone.0185056

<a name="bbmergesf">[2]</a> https://sourceforge.net/projects/bbmap/

<a name="bfcounter">[3]</a> Melsted, P. and Pritchard, J.K.: Efficient counting of k-mers in DNA sequences using a bloom filter.BMC Bioinformatics 2011 12:333.

<a name="bfcountergithub">[4]</a> https://github.com/pmelsted/BFCounter

<a name="igsimulator">[5]</a> Yana Safonova, Alla Lapidus, Jennie Lill, IgSimulator: a versatile immunosequencing simulator, Bioinformatics, Volume 31, Issue 19, 1 October 2015, Pages 3213–3215, https://doi.org/10.1093/bioinformatics/btv326

<a name="igsimulatorgithub">[6]</a> https://github.com/yana-safonova/ig_simulator

<a name="art">[7]</a> Huang W, Li L, Myers JR, Marth GT. ART: a next-generation sequencing read simulator. Bioinformatics. 2012;28(4):593–594. doi:10.1093/bioinformatics/btr708

<a name="artdl">[8]</a> https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm

<a name="seqtk">[9]</a> https://github.com/lh3/seqtk