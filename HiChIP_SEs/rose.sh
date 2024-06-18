#!/usr/bin/env bash
set PYTHONPATH /home/ben/repos/ROSE/lib
export PYTHONPATH

#awk -F \\t '{print $1 "\t" NR "\t\t" $2 "\t" $3 "\t\t.\t\t" NR}' macs_input/mef_h3k27ac_peaks.narrowPeak > rose_enhancers.gff &&


#ROSE_main.py -i rose_enhancers.gff -r ~/bio/mef_se/macs_input/BAMS/SRR17624596.bam -o rose_out --custom /home/ben/repos/ROSE/annotation/mm10_refseq.ucsc -t 2500

ROSE_main.py -i enhancer_tf.bed -r ~/bio/mef_se/macs_input/BAMS/SRR17624596.bam -o rose_out/tf --custom /home/ben/repos/ROSE/annotation/mm10_refseq.ucsc -t 2500 -c ~/bio/mef_se/macs_input/BAMS/SRR17624584.bam
