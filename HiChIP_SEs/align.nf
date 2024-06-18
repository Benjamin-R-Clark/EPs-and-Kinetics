#!/usr/bin/env nextflow


process align {
input:
path fq

output:
path "${fq.simpleName}.bam"

"""
export BOWTIE2_INDEXES=~/bio/genomes
#align using bowtie2
bowtie2 -x mm10 -p 10 -U $fq | samtools view -b -h -F 3844 -q 10 | samtools sort -o ${fq.simpleName}.bam
"""
}


process callpeaks {
publishDir "peak_result/"
input:
path bam
val WCE

output:
path "*.narrowPeak"

"""
macs2 callpeak -t $bam -c $WCE --keep-dup all -n ${bam.simpleName} -p 1e-5 -g 'mm' -f 'BAM'
"""
}

workflow {

def reads = Channel.fromPath("./mESC_TFs/*.fastq.gz")
def ctl_ch = Channel.value("~/bio/med1-mESC/SRR620141_sort.bam")

callpeaks(align(reads), ctl_ch)


}
