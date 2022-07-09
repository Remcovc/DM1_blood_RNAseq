#!/bin/bash
#parameters
#$ -o /log_o/
#$ -e /log_e/

dir=/yourdir/

python /yourdir/rmats.py \
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
