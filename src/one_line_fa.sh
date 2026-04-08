#!/usr/bin/env bash
awk 'BEGIN{RS=">"; ORS=""} NR>1{split($0,a,"\n"); header=a[1]; seq=substr($0, length(header)+2); gsub(/\n/,"",seq); print ">"header"\n"seq"\n"}' $1 > $2