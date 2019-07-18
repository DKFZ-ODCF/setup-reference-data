#!/usr/bin/env bash

###############################################################################
# download replication timings
###############################################################################
(
	mkdir_cd databases/ENCODE

	EXPECTED_MD5SUM=2a63b34a737383af2a3f7eb32801a5fa
	check_md5sum && exit 0 || echo downloading replication timing file....

	wget -c https://raw.githubusercontent.com/eilslabs/ACEseqWorkflow/github/installation/ReplicationTime_10cellines_mean_10KB.Rda

	check_md5sum
)
