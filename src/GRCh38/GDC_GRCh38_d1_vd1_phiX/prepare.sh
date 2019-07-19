#!/bin/bash

SCRIPT_DIR=$(dirname "$BASH_SOURCE")
LIB_DIR=$(readlink -f "$SCRIPT_DIR/../../scripts")


BWA_VERSION="0.7.15"

source activate setup-reference-data

#module load "bwa/$BWA_VERSION"
#module load "samtools/1.2"

set -uvex
set -o pipefail

source "$LIB_DIR/lib.sh"

# This function is used to create reduced stats files containing only the most important chromosomes,
# excluding e.g. the mitochrondrial genome or unplaced or alt-scaffold chromosomes. Note that it actually
# is executed on a tabular input containing the chr-stripped chromosome names as first column.
grepRealChromosomes() {
    grep -P '^(\d+|X|Y)\s'
}

download_GDC_GRCh38_d1_vd1_phiX() {
    local targetFile="${1:?No target filen name}"
    local md5fileDir="${2:?No MD5 file directory given}"
    local tmp=$(mktemp -d "GDC_GRCh38_d1_vd1-XXXXX")
    cachedWget \
        "https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834" \
        "$tmp/GRCh38.d1.vd1.fa.tar.gz" \
        "$md5fileDir"
    checkMD5 "$tmp/GRCh38.d1.vd1.fa.tar.gz" "$md5fileDir/GRCh38.d1.vd1.fa.tar.gz.md5"
    tar -C "$tmp" -xf "$tmp/GRCh38.d1.vd1.fa.tar.gz"
    checkMD5 "$tmp/GRCh38.d1.vd1.fa" "$md5fileDir/GRCh38.d1.vd1.fa.md5"

    download_and_add_phiX "$tmp/GRCh38.d1.vd1.fa" "$tmp/GRCh38.d1.vd1.phiX.fa" "$md5fileDir"

    cat "$tmp/GRCh38.d1.vd1.phiX.fa" > "$targetFile"  ## cat for the case "$targetFile" is /dev/stdout.
    rm -rf "$tmp"
}


initialize

primaryAssemblyName="GRCh38"
referenceGenomeName="GRCh38_d1_vd1_phiX"
outputDir="$primaryAssemblyName/$referenceGenomeName"

mkdir -p "$outputDir/$referenceGenomeName"
download_GDC_GRCh38_d1_vd1_phiX /dev/stdout "$SCRIPT_DIR" \
    | removeChrPrefixes \
    > "$outputDir/$referenceGenomeName.fa"

md5sum $(readlink -f "$outputDir/$referenceGenomeName.fa") \
    > "$SCRIPT_DIR/$referenceGenomeName.fa.md5"

bwaIndex "$outputDir/$referenceGenomeName.fa" "$outputDir/bwa-$BWA_VERSION/$referenceGenomeName.fa"
fastaIndex "$outputDir/$referenceGenomeName.fa"
statsFiles "$outputDir/$referenceGenomeName.fa" "$outputDir/stats" "grepRealChromosomes"

cleanUp