#!/bin/bash
#parameters
#$ -o /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log
#$ -e /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log
#$ -q all.q@noggo.umcn.nl


dir=/mnt/xomics/remcovc/mapping_quantification_recognition/for_paper

#/opt/cmbi/bin/STAR \
#--runThreadN 32 \
#--runMode genomeGenerate \
#--genomeDir $dir/genome_files/genomeDirSTAR_2.0.7f \
#--genomeFastaFiles $dir/genome_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#--sjdbGTFfile $dir/genome_files/Homo_sapiens.GRCh38.95.gtf \
#--sjdbOverhang 149

for sample in \
051_TGATACGC-ACTCGTTG_L003 001_GATACTGG-GAAGGTTC_L003 030_TCAGGCTT-CGCATGAT_L004 050_ATTGCGTG-GCCATAAC_L003 002_ATTCGAGG-GAGCTTGT_L003 005_CTTGTCGA-GAACATCG_L003 017_GATTACCG-CAGTCCAA_L004 016_GTCGAAGA-CAATGTGG_L004 046_GTGAGCTT-TTGATCCG_L004

do
qsub $dir/ReCognitION_complete_analysis.sh "$sample"
done