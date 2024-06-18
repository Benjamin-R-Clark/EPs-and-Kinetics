#!/bin/bash 

awk -F'\t' '{print $1, $4, " ", $2, $3, " " , ".", " ", $4}' enhancers.bed > mef_cont_enh.bed
