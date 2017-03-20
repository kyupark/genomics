#!/bin/bash -x

if [ $# -ne 3 ]
then
    echo "Usage - $0 <WGS format file> <cTea format file> <range from midpoint>"
    exit 1
fi

# Mathced TEAQ and corresponding TEA
cat $2 | xargs -n 1 -d '\n' -I % awk -v teaqline='%' 'BEGIN{FS="\t"; split(teaqline,tq); chr="chr" tq[1] } { mid=($7+$8)/2; if(chr == $2 && tq[2]+tq[6] >= mid-$3 && tq[2] <= mid+$3 ) { print teaqline; print $0 }}' $2

# Every TEAQ and matched TEA
# cat $TEAQ | xargs -n 1 -d '\n' -I % awk -v teaqline="%" 'BEGIN{FS="\t"; print teaqline; split(teaqline,tq); chr="chr" tq[1] } { if(chr == $2 && tq[2] > $3 && tq[2] < $4 ) print $0}' $TEA


#FAMILY="Alu L1 SVA PolyA"

# awk '{output= FILENAME "." $10; print >> output; close(output)}' $TEA
# awk '{output= FILENAME "." $5; print >> output; close(output)}' $TEAQ


# Matched TEAQ and TEA by Family
# for F in $FAMILY
# do
#   cat $TEAQ.$F* | xargs -n 1 -d '\n' -I % awk -v teaqline="%" 'BEGIN{FS="\t"; split(teaqline,tq); chr="chr" tq[1] } { if(chr == $2 && tq[2] > $3 && tq[2] < $4 ) { print teaqline; print $0 }}' $TEA.$F > compare-by-family.$F
# done


