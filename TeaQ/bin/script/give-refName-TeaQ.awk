#!/bin/awk -F
BEGIN { OFS="\t" }
NR == FNR { rep[$1] = $2; next }
{ $1=rep[$1] } 1

