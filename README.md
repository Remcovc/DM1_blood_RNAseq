This folder contains scripts that were used to generate the data in the sub-folders. 
 
ReCognitION_complete_analysis.sh
  This script contains the analysis from paired fastq files to count tables. It requires an input sample (base file) name as $1, which can be fed using start_ReCognitION_complete_analysis.sh. This script subsequently executes trim_galore, STAR, samtools index, umi-tools dedup, samtools index and ht-seq count. See the script file and output folders for more details
  
  rMATS.sh performs splicing analysis based on bams generated in ReCognitION_complete_analysis.sh
  
  infer_experiment.sh infers the strandedness based via RSeQC 
  
  The folder R_scripts contains R scripts of statistical modelling and to create figures
  
  The folder bams contains statistics from mapping and deduplication of RNA-seq data generated in ReCognitION_complete_analysis.sh
   
  The folder counts contains gene counts generated generated with HTSeq in ReCognitION_complete_analysis.sh 
  
  The folder trimmed_reads contains statistics from the trimming of the fastq-files using cutadapt in ReCognitION_complete_analysis.sh

  
  Raw data is available via EGA (www.ega-archive.org) under accession number .....
