This folder contains scripts to perform the raw data analysis that was used as basis for statistical modelling.


infer_experiment.sh 
    
    was used to check the strandedness of the RNA-seq data

rmats.sh  
    
    was used generate PSI counts that are analysed further in R_scripts/Figures/FigS4_PSI_CTG.R

trim_map_dedup_count.sh 
   
   was used to trim adapters etc. from the raw fastq-data, after which the reads were mapped, deduplicated and counted per gene
