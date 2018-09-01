#!/usr/bin/env bash


calling_dir=$PWD


while [ $# -ne 0 ] ; do
	cd $calling_dir
	echo "Processing :${1}:"

#	mkdir -p grouped/${1}
#	cd grouped/${1}
	mkdir -p ${1}.grouped
	cd ${1}.grouped

	cmd="seurat_group.R $1"
	echo $cmd
	$cmd

	shift
done
