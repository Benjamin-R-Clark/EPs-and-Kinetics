#!/usr/bin/env bash




macs2 callpeak -t SRR17624596.bam SRR17624597.bam \
SRR17624598.bam SRR17624599.bam SRR17624600.bam SRR17624601.bam \
-c SRR17624584.bam SRR17624585.bam \
-g 'mm' -n "mef_h3k27ac" --nomodel -p  1e-9 \
--outdir ~/bio/mef_se/macs_input/
