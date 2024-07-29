## Smartseq3 Template Workflow
### Index mm10 genome

```bash
#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --job-name=STAR_index


module load star/2.7.9a

STAR --runThreadN 4 \
 --runMode genomeGenerate \
 --genomeDir /home/clarkb/scratch/burst/ss3_2/mm10  \
 --genomeFastaFiles /home/clarkb/scratch/burst/ss3_2/mm10/mm10.fa \
 --sjdbGTFfile /home/clarkb/scratch/burst/ss3_2/mm10/mm10.ensGene.gtf  \
 --sjdbOverhang 100
```

### Create N-masked BL6/CAST genome

```bash
#!/bin/bash
#SBATCH --job-name=n_mask
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --time=00:30:00

perl /home/clarkb/projects/def-robertf/clarkb/SNPsplit/SNPsplit_genome_preparation \
 --vcf_file ~/projects/def-robertf/clarkb/Smart-seq3/allele_level_expression/CAST.SNPs.validated.vcf.gz \
 --reference_genome /home/clarkb/scratch/burst/ss3_2/mm10 --strain CAST_EiJ
```

Merge chromosome files into one:
~/scratch/burst/ss3_2/CAST/CAST_EiJ_N-masked
```
cat *.fa >> CAST_N-masked.fasta
```

### Index N-masked genome

```bash
#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --job-name=STAR_index


module load star/2.7.9a

STAR --runThreadN 4 \
 --runMode genomeGenerate \
 --genomeDir /home/clarkb/scratch/burst/ss3_2/CAST/index  \
 --genomeFastaFiles /home/clarkb/scratch/burst/ss3_2/CAST/CAST_EiJ_N-masked/CAST_N-masked.fasta \
 --sjdbGTFfile /home/clarkb/scratch/burst/ss3_2/mm10/mm10.ensGene.gtf  \
 --sjdbOverhang 100
```

### zUMIs

bash file: 


```bash
#!/bin/bash
#SBATCH --job-name=zUMI_run_2021_2h1_%j
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mail-user=benjamin_r.clark@live.com
#SBATCH --mail-type=END,FAIL

module load StdEnv/2020
module load star
module load gcc/9.3.0
module load hdf5/1.12.1
module load r/4.0.2
module load python
source ~/projects/def-robertf/clarkb/pysam-env/env/bin/activate

/home/clarkb/projects/def-robertf/clarkb/zUMIs/zUMIs.sh -y /home/clarkb/scratch/ss3_2/scripts/zUMIs_2021/zUMIs_master2.yaml
```

yaml parameters:

```
project: zUMIs_2021_star2.5.4b_FULL_h1
sequence_files:
  file1:
    name: /home/clarkb/scratch/ss3_2/reads/plate1.read1.fq.gz
    base_definition:
    - cDNA(23-100)
    - UMI(12-19)
    find_pattern: ATTGCGCAATG
  file2:
    name: /home/clarkb/scratch/ss3_2/reads/plate1.read2.fq.gz
    base_definition:
    - cDNA(1-100)
    - BC(103-108,111-118)
reference:
  STAR_index: /home/clarkb/scratch/ss3_2/CAST/N_mask_index_no_gtf
  GTF_file: /home/clarkb/scratch/ss3_2/mm10/mm10.ensGene.gtf
  additional_STAR_params: --limitSjdbInsertNsj 2000000  --clip3pAdapterSeq
    CTGTCTCTTATACACATCT
  additional_files: ~
out_dir: /home/clarkb/scratch/ss3_2/zUMIs/fullh1
num_threads: 30
mem_limit: 50
filter_cutoffs:
  BC_filter:
    num_bases: 1
    phred: 20
  UMI_filter:
    num_bases: 1
    phred: 20
barcodes:
  barcode_num: null
  barcode_file: /home/clarkb/scratch/ss3_2/reads/barcodes.txt
  automatic: no
  BarcodeBinning: 0
  nReadsperCell: 100
counting_opts:
  introns: yes
  downsampling: '0'
  strand: 0
  Ham_Dist: 1
  velocyto: no
  primaryHit: yes
  twoPass: no
make_stats: yes
which_Stage: Filtering
Rscript_exec: Rscript
STAR_exec: /home/clarkb/projects/def-robertf/clarkb/STAR-2.7.3a/source/STAR
pigz_exec: pigz
samtools_exec: /home/clarkb/projects/def-robertf/clarkb/samtools/bin/samtools
zUMIs_directory: /lustre06/project/6002165/clarkb/zUMIs
read_layout: PE
```

Notes:

The main takeaway here is that the memory limit parameter in the yaml doesn't seem to do much. When testing mem consumnption it actually takes double the amount (e.i here I gave it 50, which means it takes over 100G). Make sure the slurm job has more than double the memory requirement. This job took about 16hrs so 50G is fine.


### Merging zUMIs Output
In case you have multiple plates and want to merge them. Putting them together is not that intuitive.

First in R we'll read the UMI counts that are stored in sparse matrices. 

```R
plate1 <- readRDS("zUMIs_2021_star2.5.4b_FULL_h1.dgecounts.rds")
plate2 <- readRDS("zUMIs_2021_star2.5.4b_FULL2_h1.dgecounts.rds")

merge.counts <- function(plate1_counts, plate2_counts) {
  plate1_counts <- as.data.frame(as.matrix(plate1_counts))
  plate2_counts <- as.data.frame(as.matrix(plate2_counts))
  merged.plates <- as.matrix(merge(plate1_counts, plate2_counts, by = 0, all = TRUE))
  rnames <- merged.plates[,1]
  rownames(merged.plates) <- rnames
  merged.plates <- merged.plates[,-1]
  class(merged.plates) <- "numeric"


  out <- as(merged.plates, "dgCMatrix")

  return(out)
}

recursive.merge <- function (plate1, plate2, f){
  if(is.list(plate1) & is.list(plate2)){
    return(Map(recursive.merge, plate1, plate2, MoreArgs = list(f)))
  }
  else {
    return(f(plate1, plate2))
  }
}

merged.RDS <- recursive.merge(plate1, plate2, f = merge.counts)

```

Next is the barcodes files. Kept barcodes with UMI counts are outputted in csv. Heres a sample of how I merged them. 

```bash
(head -n 1 ../fullh1/zUMIs_output/zUMIs_2021_star2.5.4b_FULL_h1kept_barcodes.txt&& tail -n +2 ../fullh1/zUMIs_output/zUMIs_2021_star2.5.4b_FULL_h1kept_barcodes.txt && tail -n +2  ../fullh2/zUMIs_output/zUMIs_2021_star2.5.4b_FULL2_h1kept_barcodes.txt) > combined.csv
```

Re-index

```bash
awk -F , 'NR==1 {print $0; next} NR>1 {print $1, $2, NR-1}' combined.csv  > merged_barcodes.txt
```

Don't forget to merge and index the bam files.

```bash
#!/bin/bash
#SBATCH --job-name=merge_zUMIs      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=benjamin_r.clark@live.com    # Where to send mail
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=10           # Number of CPU cores per task
#SBATCH --mem=2gb                    # Job memory request
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=merge_%j.log     # Standard output and error log
module load samtools/1.16.1

samtools merge merged.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam \
../fullh1/zUMIs_2021_star2.5.4b_FULL_h1.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam \
../fullh2/zUMIs_2021_star2.5.4b_FULL2_h1.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam \
--threads 10 -f | samtools index
```


### Allele Assignment

```bash
#!/bin/bash
#SBATCH --job-name=allele_assign      # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=benjamin_r.clark@live.com           # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=24gb                    # Job memory request
#SBATCH --time=03:00:00               # Time limit hrs:min:sec
#SBATCH --output=allele_assign_j.log # Standard output and error log
#SBATCH --cpus-per-task=1             # Multi-threading

module load r/4.1.2

Rscript ~/projects/def-robertf/clarkb/Smart-seq3/allele_level_expression/get_variant_overlap_CAST.R --yaml ~/scratch/burst/ss3_2/scripts/zUMIs_master_2.run.yaml \
        --vcf  ~/projects/def-robertf/clarkb/Smart-seq3/allele_level_expression/CAST.SNPs.validated.vcf.gz
```

Notes:

Requires to use the .run.yaml file for the extra line denoting read type.


### Kinetics Estimation

```bash
#SBATCH --mail-user=benjamin_r.clark@live.com           # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=24gb                     # Job memory request
#SBATCH --time=07:00:00               # Time limit hrs:min:sec
#SBATCH --output=txburst_pipe_%j.log   # Standard output and error log
#SBATCH --cpus-per-task=30             # Multi-threading

module load scipy-stack
module load nextflow
source /home/clarkb/projects/def-robertf/clarkb/txburst/env/bin/activate
resume=false

while getopts 'r' flag
do
    case "${flag}" in
        r) resume=true ;;
        *) echo 'Unknown flag' >&2
           exit 1
    esac
done

if "$resume"; then
        echo 'resuming...'
        nextflow run main.nf -resume
else
        nextflow run main.nf
fi
```




