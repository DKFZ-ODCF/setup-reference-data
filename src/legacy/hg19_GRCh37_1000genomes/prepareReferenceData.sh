#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

outputDirectory="${1:?No output directory given}"




STAR --runMode genomeGenerate --runThreadN 8 --genomeDir "$outputDirectory" --genomeFastaFiles "$outputDirectory/GRCm38mm10_PhiX_hD3A.fa" --sjdbGTFfile hd3a_gtf/gencode-v27_hd3a.gtf --sjdbOverhang $(($readLength - 1))



==> /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation_plain.1KGRef.fa.README <==
gffread -w 1KGRef_Gencode19.fa -g hs37d5.fa gencode.v19.annotation_plain.gtf


==> /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation_plain.1KGRef.gc.README <==
perl /homes/ishaque/scripts/get_gc_content.pl gencode.v19.annotation_plain.1KGRef.fa > gencode.v19.annotation_plain.1KGRef.gc

==> /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation_plain.dexseq.gff.README <==
wget https://github.com/vivekbhr/Subread_to_DEXSeq/archive/master.zip
mv master.zip Subread_to_DEXSeq.zip
unzip Subread_to_DEXSeq.zip
python Subread_to_DEXSeq-master/dexseq_prepare_annotation2.py -f gencode.v19.annotation_plain.dexseq.gtf gencode.v19.annotation_plain.gtf gencode.v19.annotation_plain.dexseq.gff

==> /icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation_plain.transcripts.autosomal_transcriptTypeProteinCoding_nonPseudo.gtf.README <==
grep "transcript_type \"protein_coding\"" gencode.v19.annotation_plain.transcripts.gtf | grep -v ^X | grep -v ^Y | grep -v ^MT | grep -v pseudogene > gencode.v19.annotation_plain.transcripts.autosomal_transcriptTypeProteinCoding_nonPseudo.gtf
