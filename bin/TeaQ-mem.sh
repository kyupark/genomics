#!/bin/bash -x

module load seq/bwa/0.7.8

bam=$1
ref=$2
debugfolder="$1-tea-debug"
f="$debugfolder/04-contigs-f"
r="$debugfolder/04-contigs-r"

bwa mem "$ref" "$f" > "$f-mem"
bwa mem "$ref" "$r" > "$r-mem"

/home/el114/kyu/bin/script/teaify-0.1.awk "$f-mem" "$r-mem" > "$debugfolder/05-combined-mem"


