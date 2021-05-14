#!/bin/bash
#parameters
#$ -o /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log
#$ -e /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log
#$ -q all.q@narrativum.umcn.nl

# define the input parameter as "sample"
sample=$1

# define the base directory
dir=/mnt/xomics/remcovc/mapping_quantification_recognition/for_paper

# define the number of cores to use for all tools that support multithreading
cores=1

# Some notifications 
 echo "the script started with sample:"
 echo "${sample}"
 echo "at $(date +'%Y%m%d')_$(date +"%T")"

# Trim adapters and low-quality read using Trim galore on default settings. 
/mnt/home2/remcovc/bin/TrimGalore-0.6.6/trim_galore \
--paired \
--fastqc \
--gzip \
--cores $cores \
--path_to_cutadapt /home/remcovc/.local/bin/cutadapt \
-o $dir/trimmed_reads/ \
/mnt/xomics/remcovc/raw_data_recognition_2020-03/umi_extracted/"$sample"_R1.processed.fastq.gz \
/mnt/xomics/remcovc/raw_data_recognition_2020-03/umi_extracted/"$sample"_R3.processed.fastq.gz

# Map using STAR on default settings, while generating coordinate-sorted output to facilitate UMI-based deduplication and gene counting
/opt/cmbi/bin/STAR \
--runThreadN $cores \
--genomeDir /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/genome_files/genomeDirSTAR_2.7.0f/ \
--readFilesIn $dir/trimmed_reads/"$sample"_R1.processed_val_1.fq.gz $dir/trimmed_reads/"$sample"_R3.processed_val_2.fq.gz \
--readFilesCommand gunzip -c \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $dir/bams/"$sample" 

# Index the bam file for use in umi-tools
/opt/cmbi/bin/samtools index -@ $cores $dir/bams/"$sample"Aligned.sortedByCoord.out.bam 

# Deduplicate based on UMIs. Used parameters dictate the elaborate form of statistics to be output and the --spliced-is-unique parameter indicates that "two reads that start in the same position on the same strand and having the same UMI to be considered unique if one is spliced and the other is not."
/home/remcovc/.local/bin/umi_tools dedup \
-I $dir/bams/"$sample"Aligned.sortedByCoord.out.bam \
--spliced-is-unique \
--paired \
-L $dir/bams/dedup/"$sample".dedup_logfile.txt \
--output-stats=$dir/bams/dedup/"$sample".dedup_stats \
-S $dir/bams/dedup/"$sample".dedup.bam

# Index the deduplicated bam file for htseq
/opt/cmbi/bin/samtools index -@ $cores $dir/bams/dedup/"$sample".dedup.bam

# Gene counting using HTSeq, default settings, same genome version/gtf file used as in STAR. Stranded=yes because rseqc infer_experiment explained most reads (85-95%) with "1++,1--,2+-,2-+". See the infer_experiment.sh script and the output in bams/dedup/infer_experiment
python -m HTSeq.scripts.count \
$dir/bams/dedup/"$sample".dedup.bam \
/mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/genome_files/Homo_sapiens.GRCh38.95.gtf \
--format bam \
--order pos \
--stranded=yes \
> $dir/counts/"$sample".counts.txt 