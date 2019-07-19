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

# To remove HLA loci and alt-scaffold locus sequences from this assembly file.
grepNonAltChromosomes() {
    "$LIB_DIR/grepChromosomes.py" --nomatch '(?:_alt$|^HLA)'
}

# This function is used to create reduced stats files containing only the most important chromosomes,
# excluding e.g. the mitochrondrial genome or unplaced or alt-scaffold chromosomes. Note that it actually
# is executed on a tabular input containing the chr-stripped chromosome names as first column.
grepRealChromosomes() {
    grep -P '^(\d+|X|Y)\s'
}


initialize

primaryAssemblyName="GRCh38"
referenceGenomeName="GRCh38_decoy_ebv_phiX"
outputDir="$primaryAssemblyName/$referenceGenomeName"

mkdir -p "$outputDir/$referenceGenomeName"
download_core_ref_GRCh38_hla_decoy_ebv_phiX "$outputDir/${referenceGenomeName}_alt_hla.fa" "$SCRIPT_DIR"

cat "$outputDir/${referenceGenomeName}_alt_hla.fa" \
    | grepNonAltChromosomes \
    | removeChrPrefixes \
    > "$outputDir/$referenceGenomeName.fa"
rm "$outputDir/${referenceGenomeName}_alt_hla.fa"

md5sum $(readlink -f "$outputDir/$referenceGenomeName.fa") \
    > "$SCRIPT_DIR/$referenceGenomeName.fa.md5"

bwaIndex "$outputDir/$referenceGenomeName.fa" "$outputDir/bwa-$BWA_VERSION/$referenceGenomeName.fa"
fastaIndex "$outputDir/$referenceGenomeName.fa"
statsFiles "$outputDir/$referenceGenomeName.fa" "$outputDir/stats" "grepRealChromosomes"

cleanUp
