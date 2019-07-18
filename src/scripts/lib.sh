#!/bin/bash

initialize() {
    CACHE_DIR="$PWD/cache"
    mkdir -p "$CACHE_DIR"
}


cleanUp() {
    echo "Cache directory ($CACHE_DIR) needs to be cleaned. Consider also removing old temp-directories." >> /dev/stderr
}


cachedWget() {
    local url="${1:?No URL given}"
    local targetFile="${2:?No target filename given}"
    local md5fileDir="${3:?No MD5 file dir given}"
    local cacheDir="${4:-${CACHE_DIR:?No cache directory}}"
    local cachedFile="$cacheDir/"$(basename "$targetFile")
    if [[ ! -f "$cachedFile" ]]; then
        wget -c "$url" -O "$cachedFile"
    else
        echo "Cached file '$cachedFile' exists. Skipping re-download ..." >> /dev/stderr
    fi
    cp "$cachedFile" "$targetFile"
}

# This requires that the MD5 file contains a line with the MD5 sum and a tab-separated '-', representing the standard error.
checkMD5() {
    local file="${1:?No file to check given}"
    local md5file="${2:?No MD5 file given}"
    cat "$file" | md5sum --strict -c "$md5file" >> /dev/null
}


download_phiX() {
    local targetFile="${1:?No target file name}"
    local md5fileDir="${2:?No MD5 file directory given}"
    local tmp=$(mktemp -d "phiX-XXXXX")
    cachedWget \
        "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz" \
        "$tmp/PhiX_Illumina_RTA.tar.gz" \
        "$md5fileDir"
    checkMD5 "$tmp/PhiX_Illumina_RTA.tar.gz" "$md5fileDir/PhiX_Illumina_RTA.tar.gz.md5"
    tar -C "$tmp" -xf "$tmp/PhiX_Illumina_RTA.tar.gz" "PhiX/Illumina/RTA/Sequence/Chromosomes/phix.fa"
    cat "$tmp/PhiX/Illumina/RTA/Sequence/Chromosomes/phix.fa" > "$targetFile"  ## cat for the case "$targetFile" is /dev/stdout.
    checkMD5 "$targetFile" "$md5fileDir/PhiX_Illumina_RTA.fa.md5"
    rm -rf "$tmp"
}

download_core_ref_GRCh38_hla_decoy_ebv() {
    local targetFile="${1:?No target filen name}"
    local md5fileDir="${2:?No MD5 file directory given}"
    local tmp=$(mktemp -d "core_ref_GRCh38_hla_decoy_ebv-XXXXX")
    cachedWget \
        "ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz" \
        "$tmp/core_ref_GRCh38_hla_decoy_ebv.tar.gz" \
        "$md5fileDir"
    checkMD5 "$tmp/core_ref_GRCh38_hla_decoy_ebv.tar.gz" "$md5fileDir/core_ref_GRCh38_hla_decoy_ebv.tar.gz.md5"
    tar -C "$tmp" -xf "$tmp/core_ref_GRCh38_hla_decoy_ebv.tar.gz" "core_ref_GRCh38_hla_decoy_ebv/genome.fa"
    cat "$tmp/core_ref_GRCh38_hla_decoy_ebv/genome.fa" > "$targetFile"  ## cat for the case "$targetFile" is /dev/stdout.
    checkMD5 "$targetFile" "$md5fileDir/core_ref_GRCh38_hla_decoy_ebv.fa.md5"
    rm -rf "$tmp"
}

download_and_add_phiX() {
    local sourceFile="${1:?No source file name given}"
    local targetFile="${2:?No target file name given}"
    local md5fileDir="${3:?No MD5 file directory given}"
    local tmp=$(mktemp -d "fasta-XXXXX")
    download_phiX "$tmp/phiX.fa" "$md5fileDir"
    cat "$sourceFile" "$tmp/phiX.fa" > "$targetFile"  ## cat for the case "$targetFile" is /dev/stdout.
    rm -rf "$tmp"
}

download_core_ref_GRCh38_hla_decoy_ebv_phiX() {
    local targetFile="${1:?No target file name given}"
    local md5fileDir="${2:?No MD5 file directory given}"
    local tmp=$(mktemp -d "fasta-XXXXX")
    download_core_ref_GRCh38_hla_decoy_ebv "$tmp/core_ref_GRCh38_hla_decoy_ebv.fa" "$md5fileDir"
    download_and_add_phiX "$tmp/core_ref_GRCh38_hla_decoy_ebv.fa" "$targetFile" "$md5fileDir"
    rm -rf "$tmp"
}

# The "prefix" should be a path plus an index name (which could be for instance the FASTA file name).
# The actual output path will be determined with dirname from the prefix.
bwaIndex() {
    local fasta="${1:?No FASTA file given}"
    local prefix="${2:-$fasta}"
    local indexDir=$(dirname "$prefix")
    local indexName=$(basename "$prefix")
    local absoluteFasta=$(readlink -f "$fasta")
    local fastaFileName=$(basename "$fasta")
    mkdir -p "$indexDir"
    (
        ln -sf "$absoluteFasta" "$indexDir"
        pushd "$indexDir"
        bwa index -p "$indexName" "$fastaFileName"
    )
}

fastaIndex() {
    local fasta="${1:?No FASTA file given}"
    echo "$PWD"
    samtools faidx "$fasta"
}

# Create a base count file only counting A, C, T, G and a variant only for the main chromosomes {1..22,X,Y}. The variant is only for the Roddy
# alignment workflow and constitutes a reduced set on which analyses are done.
statsFiles() {
    local fasta="${1:?No FASTA file given}"
    local outputDir="${2:?No output directory given}"
    local realChromosomesFunName="${3:?No function name to grep real chromosomes with stream filtering}"
    mkdir -p "$outputDir"
    local base=$(basename "$fasta")
    local ciFile="$outputDir/$base.contigInformation.tsv"

    "$LIB_DIR/refGenomeBaseCounts.py" "$fasta" \
        > "$ciFile"

    # Full lengths of restricted chromosome set
    cat "$ciFile" \
        | cut -f 1,2 \
        | "$realChromosomesFunName" \
        > "$outputDir/$base.chrLength.tsv"

    # Non-N lengths of all chromosomes and contigs
    local fullACGTFile="$outputDir/$base.chrLenOnlyACGT.tsv"
    cat "$ciFile" \
        | cut -f 1,3 \
        > "$fullACGTFile"

    # Non-N lengths of restricted chromosomes set
    local restrictedACGTFile="$outputDir/$base.chrLenOnlyACGT_realChromosomes.tsv"
    cat "$fullACGTFile" \
        | "$realChromosomesFunName" \
        > "$restrictedACGTFile"
}

removeChrPrefixes() {
    local infasta="${1:-/dev/stdin}"
    local outfasta="${1:-/dev/stdout}"
    cat "$infasta" \
        | sed -r 's/>chr/>/' \
        > "$outfasta"
}

## Create a .dict file that looks like a SAM header for use with GATK
#gatkDict() {
#
#    cd /icgc/ngs_share/assemblies/legacy.hg_GRCh38
#    echo -e "@HD\tVN:1.0\tSO:unsorted" > sequence/legacy.hg_GRCh38.dict
#    sed '1d' stats/legacy.hg_GRCh38.fa.contigInformation.txt | awk '{print "@SQ\tSN:" $1 "\tLN:" $5 "\tUR:file:/icgc/ngs_share/assemblies/legacy.hg_GRCh38/sequence/legacy.hg_GRCh38.fa"}' >> sequence/legacy.hg_GRCh38.dict
#    cd sequence
#    ln -s legacy.hg_GRCh38.dict Homo_sapiens.GRCh38.dna_sm.primary_assembly.dict
#    cd ../../../indexes/bwa/bwa06
#    ln -s ../../../sequence/Homo_sapiens.GRCh38.dna_sm.primary_assembly.dict legacy.hg_GRCh38.dict
#}


## Tell Bash to abort on error.
#set -o pipefail
#set -u
#set -eE
#trap 'echo "Download incomplete. Please restart script."' ERR
#
#
## convenience function to create a directory and change into it
#mkdir_cd() {
#	mkdir -p "$1"
#	cd "$1"
#}
#
## compute an MD5 sum over all files found recursively in the current directory
## and check it against the MD5 sum given in the variable EXPECTED_MD5SUM
#check_md5sum() {
#	local FILES=$(find -type f | sort)
#	local MD5SUM=$([ -n "$FILES" ] && cat $FILES | md5sum | cut -f1 -d' ')
#	[ "$EXPECTED_MD5SUM" = "$MD5SUM" ]
#}
#
#
#
#
#downloadReferenceGenome() (
#    local fastaGzUrl="${1:-ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz}"
#    local expectedMd5="${2:-12a0bed94078e2d9e8c00da793bbc84e}"
#	EXPECTED_MD5SUM=$expectedMd5
#
#	mkdir_cd sequence/1KGRef
#
#	check_md5sum && exit 0 || echo downloading reference genome....
#
#	wget -c "$fastaGzUrl"
#	gunzip hs37d5.fa.gz
#
#	check_md5sum
#)
#
## CITATION
## As a NCBI Resource: "Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K. dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;29(1):308-11."
## As a whole for a specific build (use this!) : "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. (dbSNP Build ID: 141 ). Available from: http://www.ncbi.nlm.nih.gov/SNP/"
## A single or a range of Submitted SNP (ss) or Reference SNP (rs) entries: "Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center for Biotechnology Information, National Library of Medicine. dbSNP accession:{ss1 or ss1 – ss100}, (dbSNP Build ID: 141). Available from: http://www.ncbi.nlm.nih.gov/SNP/"
#downloadDbSnpDatabase() (
#    local DBSNP_VERSION="${1:-135}"
#	local DBSNP_BASE_URL="${2:-ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF}"
#	local EXPECTED_MD5SUM="${3:-fed2a31b5a5d8fe12e072576c0c17199}"
#
#	mkdir_cd databases/dbSNP/dbSNP_$DBSNP_VERSION
#
#	check_md5sum && exit 0 || echo downloading dbSNP file....
#
#	# DOWNLOAD
#	wget -c "$DBSNP_BASE_URL/README.txt"
#	wget -c "$DBSNP_BASE_URL/00-All.vcf.gz"
#
#	# POST PROCESSING
#	# extract SNPs from dbSNP version 135 and older
#	if [[ "$DBSNP_VERSION" <= 135 ]]; then
#	    zcat 00-All.vcf.gz \
#    	    | awk '/^#/{print} /VC=SNV/{ v=$8; sub(/.*dbSNPBuildID=/, "", v); sub(/;.*/, "", v); if (v~/^[0-9]+$/ && int(v)<='$DBSNP_VERSION') print }' \
#    	    | bgzip > 00-All.SNV.vcf.gz
#    	tabix -p vcf 00-All.SNV.vcf.gz
#    else
#        echo "No code present to download dbSNP versions newer than 135."
#        return 1
#    fi
#
#	# CLEANUP
#	rm -f 00-All.vcf.gz
#
#	check_md5sum
#)
#
#
#downloadMappabilityFile() (
#	local mappabilityBigWig="${1:-http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig}"
#	local EXPECTED_MD5SUM="${2:-3d12d0a4d7afdb52cfd10f886d48b5f0}"
#
#	mkdir_cd databases/UCSC
#
#	check_md5sum && exit 0 || echo downloading mappability file....
#
#	wget -c "$mappabilityBigWig"
#
#	local bigwig=$(basename "$mappabilityBigWig")
#	local bedgraph=$(basename "$bigwig" .bigwig)_chr.bedGraph.gz
#	bigWigToBedGraph "$bigwig" /dev/stdout | bgzip > "$bedgraph"
#	tabix -p bed "$bedgraph"
#	rm -f "$bigwig"
#
#	check_md5sum
#)
#
#
#downloadReplicationTimingFile() (
#	local replicationTimingRdaFile="${1:-https://raw.githubusercontent.com/eilslabs/ACEseqWorkflow/github/installation/ReplicationTime_10cellines_mean_10KB.Rda}"
#	local EXPECTED_MD5SUM="${2:-2a63b34a737383af2a3f7eb32801a5fa}"
#
#	mkdir_cd databases/ENCODE
#
#	check_md5sum && exit 0 || echo downloading replication timing file....
#
#	wget -c "$replicationTimingRdaFile"
#
#	check_md5sum
#)
#
#downloadChromosomeStatistics() (
#    local chrInfoGz="${1:-http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz}"
#    local chrInfoTxt="${2:-https://raw.githubusercontent.com/eilslabs/ACEseqWorkflow/github/installation/hg19_GRch37_100genomes_gc_content_10kb.txt}"
#	local EXPECTED_MD5SUM="${3:-801bdaa8c3b0d5c18a0637b0b29fd337}"
#
#	mkdir_cd stats
#
#	check_md5sum && exit 0 || echo downloading stats files....
#
#	wget -c "$chrInfoGz"
#
#	local infoFile=$(basename "$chrInfoGz")
#
#	zcat "$infoFile" | grep -Pv "(_)|(chrM)" | sed -e '1i\#chrom\tsize\tfileName' > chrlengths.txt
#	rm -f "$infoFile"
#
#	wget -c "$chrInfoTxt"
#
#	check_md5sum
#)
#
#downloadImputeDatabase() (
#    local allPhase1ShapeItNoMonoFile="${1:-https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz}"
#    local all1000gPhase1ImputeFile="${2:-https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz}"
#	local EXPECTED_MD5SUM="${3:-261a28d6b6917340cd82ada2d7185e17}"
#
#	mkdir_cd databases/1000genomes/IMPUTE
#
#	check_md5sum && exit 0 || echo downloading impute files....
#
#	wget -c "$allPhase1ShapeItNoMonoFile"
#    local allPhase1ShapeItNoMonoFile_bn=$(basename "$allPhase1ShapeItNoMonoFile")
#	tar -xzvf "$allPhase1ShapeItNoMonoFile_bn"
#	rm -f "$allPhase1ShapeItNoMonoFile_bn"
#
#	wget -c "$all1000gPhase1ImputeFile"
#	local all1000gPhase1ImputeFile_bn=$(basename "$all1000gPhase1ImputeFile")
#	tar -xzvf "$all1000gPhase1ImputeFile_bn"
#	rm -f "$all1000gPhase1ImputeFile_bn"
#
#	check_md5sum
#)
#
#
### Needs tabix
### Reconstructing the ancestral germ line methylation state of young repeats.
### Feuerbach L1, Lyngsø RB, Lengauer T, Hein J. Mol Biol Evol. 2011 Jun;28(6):1777-84. doi: 10.1093/molbev/msr001. Epub 2011 Jan 6.
#downloadCgiMountains() {
#    local cgiMountainsUrl="${1:-https://cgihunter.bioinf.mpi-inf.mpg.de/newdownload/CgiMountain_hg19_liftover.bed.gz}"
#    local outfile="CgiMountains_chr.bed.gz"
#    wget "$cgiMountainsUrl"
#    local bName=$(basename "$cgiMountainsUrl")
#    zcat "$bName" | tail -n +2 | cut -f 1-6 | sort -k 1d,1 -k 2n,3 | bgzip -c > "$outfile"
#    tabix $outfile
#}
#
#
#
#download
#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/ftp/technical/working/20120213_phase1_integrated_release_version1/ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf.gz
