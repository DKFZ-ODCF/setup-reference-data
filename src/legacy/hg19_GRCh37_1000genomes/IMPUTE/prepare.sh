#!/usr/bin/env bash


###############################################################################
# download IMPUTE database
###############################################################################
(
	mkdir_cd databases/1000genomes/IMPUTE

	EXPECTED_MD5SUM=261a28d6b6917340cd82ada2d7185e17
	check_md5sum && exit 0 || echo downloading impute files....

	wget -c https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz
	tar -xzvf ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz
	rm -f ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz

	wget -c https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
	tar -xzvf ALL_1000G_phase1integrated_v3_impute.tgz
	rm -f ALL_1000G_phase1integrated_v3_impute.tgz

	check_md5sum
)
