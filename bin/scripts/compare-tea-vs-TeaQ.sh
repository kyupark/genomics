#!/bin/bash -x


awk '{ OFS="\t"; print $1, $2, $3}' $1 | xargs -n 1 -d '\n' -I % sh -c 'grep "%" $1 ; grep "%" $2 ; echo ""'
