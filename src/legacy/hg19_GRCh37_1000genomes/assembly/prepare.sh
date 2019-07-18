#!/bin/bash

# tell bash to abort on error
set -o pipefail
set -u
set -eE
trap 'echo "Download incomplete. Please restart script."' ERR

# convenience function to create a directory and change into it
mkdir_cd() {
	mkdir -p "$1"
	cd "$1"
}

# compute an MD5 sum over all files found recursively in the current directory
# and check it against the MD5 sum given in the variable EXPECTED_MD5SUM
check_md5sum() {
	local FILES=$(find -type f | sort)
	local MD5SUM=$([ -n "$FILES" ] && cat $FILES | md5sum | cut -f1 -d' ')
	[ "$EXPECTED_MD5SUM" = "$MD5SUM" ]
}


###############################################################################
# create download directory
###############################################################################

mkdir_cd hg19_GRCh37_1000genomes

###############################################################################
# download reference genome
###############################################################################
(
	mkdir_cd sequence/1KGRef

	EXPECTED_MD5SUM=12a0bed94078e2d9e8c00da793bbc84e
	check_md5sum && exit 0 || echo downloading reference genome....

	wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	gunzip hs37d5.fa.gz

	check_md5sum
)


###############################################################################
# download chromosome statistics
###############################################################################
(
	mkdir_cd stats

	EXPECTED_MD5SUM=801bdaa8c3b0d5c18a0637b0b29fd337
	check_md5sum && exit 0 || echo downloading stats files....

	wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz
	zcat chromInfo.txt.gz | grep -Pv "(_)|(chrM)" | sed -e '1i\#chrom\tsize\tfileName' > chrlengths.txt
	rm -f chromInfo.txt.gz

	wget -c https://raw.githubusercontent.com/eilslabs/ACEseqWorkflow/github/installation/hg19_GRch37_100genomes_gc_content_10kb.txt

	check_md5sum
)
