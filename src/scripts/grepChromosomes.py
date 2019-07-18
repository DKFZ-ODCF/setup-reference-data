#!/usr/bin/env python3
# @author: Philip R. Kensche
# Read in a FASTA file and one or multiple patterns for chromosomes to keep and return only those with a matching entry ID.

import sys
import re
from Bio import SeqIO
import argparse

def asCompiledRegularExpressions(expressions):
    return list(map(lambda e: re.compile(e), expressions))


def getOptions(options):
    parser = argparse.ArgumentParser(description="Grep sequences from a file based on regular expressions matched to the entry identifiers.")
    parser.add_argument(
        '--match', dest="shouldMatch", type=str, action="append", default=['.*'],
        help="Regular expressions of which at least one should match."
    )
    parser.add_argument(
        "--nomatch", dest="shouldNotMatch", type=str, action="append", default=[],
        help="Regular expressions of which none should match."
    )
    args = parser.parse_args(options)
    args.shouldMatch = asCompiledRegularExpressions(args.shouldMatch)
    args.shouldNotMatch = asCompiledRegularExpressions(args.shouldNotMatch)
    return args


def searchRe(regular_expressions, string):
    '''Try matching any of the regular expressions in the input string. Return True only on success of at least one expression.'''
    for exp in regular_expressions:
        if re.search(exp, string):
            return True
    return False



try:
    args = getOptions(sys.argv[1:])
except Exception as e:
    print("Error in parameters: {e}", file=sys.stderr)
    exit(1)


in_fh = sys.stdin
out_fh = sys.stdout

format = "fasta"

kept = 0
dropped = 0
for record in SeqIO.parse(in_fh, format):
    if searchRe(args.shouldMatch, record.id) and not searchRe(args.shouldNotMatch, record.id):
        print("Keeping '{}' ...".format(record.id), file=sys.stderr)
        SeqIO.write(record, out_fh, format)
        kept += 1
    else:
        print("Dropping '{}' ...".format(record.id), file=sys.stderr)
        dropped += 1

out_fh.close()

print("Kept {} and dropped {} sequences.".format(kept, dropped), file=sys.stderr)