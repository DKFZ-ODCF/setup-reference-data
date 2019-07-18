#!/usr/bin/env bash


###############################################################################
# download mappability file
###############################################################################
(
	mkdir_cd databases/UCSC

	EXPECTED_MD5SUM=3d12d0a4d7afdb52cfd10f886d48b5f0
	check_md5sum && exit 0 || echo downloading mappability file....

	wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
	bigWigToBedGraph wgEncodeCrgMapabilityAlign100mer.bigWig /dev/stdout | bgzip > wgEncodeCrgMapabilityAlign100mer_chr.bedGraph.gz
	tabix -p bed wgEncodeCrgMapabilityAlign100mer_chr.bedGraph.gz
	rm -f wgEncodeCrgMapabilityAlign100mer.bigWig

	check_md5sum
)