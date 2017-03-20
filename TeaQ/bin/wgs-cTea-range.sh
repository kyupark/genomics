#!/bin/bash -x

if [ $# -ne 3 ]
then
    echo "Usage - $0 <WGS format file> <cTea format file> <range from midpoint>"
    exit 1
fi

# Mathced TEAQ and corresponding TEA
cat $2 | xargs -n 1 -d '\n' -I % awk -v teaqline='%' -v range=$3 'BEGIN{FS="\t"; split(teaqline,tq); chr="chr" tq[1] } { mid=($7+$8)/2; if(chr == $2 && tq[2]+length(tq[7]) >= mid-range && tq[2] <= mid+range ) { print teaqline; print $0 }}' $1


