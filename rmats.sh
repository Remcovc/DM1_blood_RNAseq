#!/bin/bash
#parameters
#$ -o /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log/
#$ -e /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log/
#$ -q all.q@narrativum.umcn.nl

dir=/mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/

python /mnt/home2/remcovc/miniconda3/envs/rMATS/rMATS/rmats.py \
--b1 $dir/bam_paths.txt  \
--gtf $dir/genome_files/Homo_sapiens.GRCh38.95.gtf \
--od $dir/rMATS/ \
--tmp $dir/rMATS/temp \
-t paired \
--readLength 150 \
--variable-read-length \
--nthread 32 \
--novelSS \
--libType fr-firststrand \
--statoff