#!/usr/bin/env python3
# @author: Philip R. Kensche
# Read in a FASTA file and return a TSV with one entry per FASTA entry and three columns (entryId, totalLength, lengthOnlyATCG).
# This is the information needed for the "stats" file.

import sys
import os.path
from Bio import SeqIO
from collections import Counter

if len(sys.argv) != 2:
    sys.exit("Please supply FASTA file, and only FASTA file! Produce a headerless TSV with columns (entryId, totalLength, lenthOnlyATCG.")
if len(sys.argv) == 2:
    if not os.path.isfile(sys.argv[1]):
        sys.exit("Supplied FASTA file cannot be found!")

in_fh = open(sys.argv[1], "r")
out_fh = sys.stdout

print("Reading sequence data from " + sys.argv[1], file=sys.stderr)

for record in SeqIO.parse(in_fh, "fasta"):
    print("Processing entry '{}' ...".format(record.id), file=sys.stderr)
    counter = Counter(record.seq)
    withN = sum(counter.values())
    withoutN = sum(list(map(lambda c: counter[c], ["A", "a", "T", "t", "C", "c", "G", "g"])))
    print("\t".join([record.id, str(withN), str(withoutN)]), file=out_fh)

in_fh.close()
