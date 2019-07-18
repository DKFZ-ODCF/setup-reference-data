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
# download dbSNP database
###############################################################################
(
	DBSNP_VERSION=135
	mkdir_cd databases/dbSNP/dbSNP_$DBSNP_VERSION

	EXPECTED_MD5SUM=fed2a31b5a5d8fe12e072576c0c17199
	check_md5sum && exit 0 || echo downloading dbSNP file....

	# CITATION
	# As a NCBI Resource: "Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K. dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;29(1):308-11."
	# As a whole for a specific build (use this!) : "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. (dbSNP Build ID: 141 ). Available from: http://www.ncbi.nlm.nih.gov/SNP/"
	# A single or a range of Submitted SNP (ss) or Reference SNP (rs) entries: "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. dbSNP accession:{ss1 or ss1 – ss100}, (dbSNP Build ID: 141). Available from: http://www.ncbi.nlm.nih.gov/SNP/"

	# DOWNLOAD
	DBSNP_BASE_URL="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF"
	wget -c "$DBSNP_BASE_URL/README.txt"
	wget -c "$DBSNP_BASE_URL/00-All.vcf.gz"

	# POST PROCESSING
	# extract SNPs from dbSNP version 135 and older
	zcat 00-All.vcf.gz |
	awk '/^#/{print} /VC=SNV/{ v=$8; sub(/.*dbSNPBuildID=/, "", v); sub(/;.*/, "", v); if (v~/^[0-9]+$/ && int(v)<='$DBSNP_VERSION') print }' |
	bgzip > 00-All.SNV.vcf.gz
	tabix -p vcf 00-All.SNV.vcf.gz

	# CLEANUP
	rm -f 00-All.vcf.gz

	check_md5sum
)


echo "All files downloaded successfully"
