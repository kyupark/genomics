#!/bin/bash -x

BAM=$1
REFNAME="$1.refName"

samtools view -H $BAM | awk 'BEGIN{i=0} { if ( $2 ~ /^SN:/ ) { gsub(/^SN:/,"",$2); print i, $2; i++ } }' > $REFNAME

