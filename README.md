This folder contains scripts that were used to generate the data in the sub-folders. 

start_ReCognitION_complete_analysis.sh
  This script is used to start the analysis in ReCognitION_complete_analysis.sh. This contains the genome generation for STAR and sample name specification
  
ReCognitION_complete_analysis.sh
  This script contains the analysis from paired fastq files to count tables. It requires an input sample (base file) name as $1, which can be fed using start_ReCognitION_complete_analysis.sh. This script subsequently executes trim_galore, STAR, samtools index, umi-tools dedup, samtools index and ht-seq count. See the script file and output folders for more details