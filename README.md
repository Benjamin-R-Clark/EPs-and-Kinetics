# EPs-and-Kinetics
Master thesis workspace focusing on the effects of enhancer-promoter contacts on transcritptional kinetics.


This space is currently under construction, I'm slowly adding md presentations and figures. 

What I'm providing here in this README is a general overview of all the tools I've used in investigating the relationship between enhancer-promoter contacts (EPs) and transcriptional kinetics. Namely, if transcription occurs in discrete, stochastic bursts, do enhancers when physically interacting with promoters change the frequency or the size of these bursts? What about super-enhancers?

To begin we need to estimate transcriptional kinetics from scRNAseq data from MEF cells. We are able to do this using public data from the Sandberg group.
### 1.) [SS3_scRNAseq_Workflow.md](https://github.com/Benjamin-R-Clark/eps-and-kinetics-thesis/blob/main/SS3_scRNAseq_Workflow.md)
  Benchling entry on scRNAseq workflow using zUMIs on a hybrid mouse strain. Mostly scripts and notes on getting things to run on an HPC, here using slurm on compute canada servers.

### 2.) [txburst_pipe](https://github.com/Clarkvale/txburst_pipe)
  Nextflow pipeline for calling transcriptional kinetics from allele-specific UMI counts. Relies on scripts from  [sandberg lab](https://github.com/sandberg-lab/txburst).

### 3.) [ChipSeq Nextflow](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/blob/main/HiChIP_SEs/align.nf)
  In order to define Super-Enhancers we first take MEF or mESC H3K27ac ChIPseq reads and align, filter and call peaks using MACS2. This is partly accomplished using this nextflow pipeline.

### 4.)[ROSE script](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/blob/main/HiChIP_SEs/rose.sh)
  We then run ROSE using the output and the control input. Constituent enhancers are defined from the individual H3K27ac peaks filtered from promoters:

### 5.)[MEF Super-Enhancer Analysis](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/blob/main/HiChIP_SEs/MEF/mef_se.md)
  We lastly need to map SEs to promoters. The final analysis is in R or R-markdown. Do Super-enhancers influence burst frequency or burst size? Or both? What about genomic distance? Do linkages actually matter?:




### Misc
[HiChIP_SEs](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/tree/main/HiChIP_SEs)
  Work repos for identfying super enhancers (SEs) from chipseq files in MEFs and mESCs. I use these annotations with H3K27ac HICHIP chromatin loops to map promoters with SEs. I then incorporate transcriptional kinetics data to determine what effect these linkages have. I have sample bash scripts and short nextflow pipelines to show how I did secondary analysis (alignment, peak calling and SE calling). 
  
[tx_noise](https://github.com/Clarkvale/txnoise/tree/74118078493861024634fedda180b60544bd8bd4)
  Workspace for calling transcriptional noise values from allele-specific UMI counts.

[hichip-burst.R](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/blob/main/hichip/hichip-burst.R)
  Script for mapping enhancers to promoters and adding the kinetics for each gene. Bootstrapped linear models are performed to see if the number of enhancers changes kinetic outputs. The short answer is yes, but in interesting ways (see poster).

[ep_networks.Rmd](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/blob/main/hichip/ep_networks.rmd)
  Markdown script for generating network graphs of EPs from H3K27ac HI-ChIP data. Initial hypothesis was to see if network interactions between enhancers and promoters using common graph statistics (such as centrality) influence kinetics of the members   within. It turns out these networks can be quite complex and difficult draw significant conclusions out of. The script does make some really pretty graph networks though. 
[enh_knn.R](https://github.com/Benjamin-R-Clark/eps-and-kinetics-thesis/blob/main/enh_knn.R)
  My attempt at implementing a random forest model on enhancer sequences to predict transcription burst frequency. Ideally the model would be able to find combinatorial TF motifs to predict burst frequency Unfortunately the model failed, most likely to the size of the training set (~150) and the size of enhancer sequeences (~100bp). A better model would isolate specific regions using ATACseq. 

  ### Project Poster

  ![BC_2024_IRCM-1](https://github.com/user-attachments/assets/9fea9f99-5c22-40ed-9239-c7647a2d4598)

  
  🚧🚧🚧🚧🚧 Under construction 🚧🚧🚧🚧🚧
