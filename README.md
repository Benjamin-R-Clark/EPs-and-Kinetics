# EPs-and-Kinetics
Master thesis workspace focusing on the effects of enhancer-promoter contacts on transcritptional kinetics.


This space is currently under construction, I'm slowly adding md presentations and figures. 

#### [HiChIP_SEs](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/tree/main/HiChIP_SEs)
  Work repos for identfying super enhancers (SEs) from chipseq files in MEFs and mESCs. I use these annotations with H3K27ac HICHIP chromatin loops to map promoters with SEs. I then incorporate transcriptional kinetics data to determine what effect these linkages have. I have sample bash scripts and short nextflow pipelines to show how I did secondary analysis (alignment, peak calling and SE calling). The final analysis is in R or R-markdown such as:

  [MEF Super-Enhancer Analysis](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/blob/main/HiChIP_SEs/MEF/mef_se.md)

#### [txburst_pipe](https://github.com/Clarkvale/txburst_pipe)
  Nextflow pipeline for calling transcriptional kinetics from allele-specific UMI counts. Relies on scripts from  [sandberg lab](https://github.com/sandberg-lab/txburst).

#### [tx_noise](https://github.com/Clarkvale/txnoise/tree/74118078493861024634fedda180b60544bd8bd4)
  Workspace for calling transcriptional noise values from allele-specific UMI counts.

#### [SS3_scRNAseq_Workflow.md](https://github.com/Benjamin-R-Clark/EPs-and-Kinetics/SS3_scRNAseq_Workflow.md)
  Benchling entry on scRNAseq workflow using zUMIs on a hybrid mouse strain. Mostly scripts and notes on getting things to run on an HPC, here using slurm on compute canada servers.


  ### Project Poster

  ![BC_2024_IRCM-1](https://github.com/user-attachments/assets/9fea9f99-5c22-40ed-9239-c7647a2d4598)

  
  🚧🚧🚧🚧🚧 Under construction 🚧🚧🚧🚧🚧
