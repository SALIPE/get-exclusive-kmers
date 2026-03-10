#!/bin/bash

PROJECTHOME=~/Desktop/get-exclusive-kmers/GEKMERS

HBV=~/Desktop/datasets/HBV/data
REF_HBV=~/Desktop/datasets/HBV/refseq.fasta

KMER=$1

cd $PROJECTHOME && julia --project src/GetExclusiveKmers.jl \
        --directory $HBV \
        --reference $REF_HBV \
        --k-len $KMER \
        --use-entropy
        # --use-gramep
