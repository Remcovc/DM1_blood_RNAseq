#!/bin/bash
#parameters
#$ -o /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log
#$ -e /mnt/xomics/remcovc/mapping_quantification_recognition/for_paper/log
#$ -q all.q@noggo.umcn.nl

# define the base directory
dir=/mnt/xomics/remcovc/mapping_quantification_recognition/for_paper

# run STAR genome generate to generate an indexed genome directory from fasta/gtf files. sjdbOverhang is at 149 here because the read lengths used were 150 bp. STAR version was 2.7.0f. Genome version GRCh38.95 from ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens and ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/
/opt/cmbi/bin/STAR \
--runThreadN 32 \
--runMode genomeGenerate \
--genomeDir $dir/genome_files/genomeDirSTAR_2.7.0f \
--genomeFastaFiles $dir/genome_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile $dir/genome_files/Homo_sapiens.GRCh38.95.gtf \
--sjdbOverhang 149

# Define sample names. For each of these names the ReCognitION_complete_analysis.sh will be started in parallel. The fastq files need to be named "sample_R1.processed.fastq.gz" and "sample_R3.processed.fastq.gz" and should contain UMIs in their read name
for sample in \
001_GATACTGG-GAAGGTTC_L003 031_CCTTGTAG-TTCCAAGG_L004 \
002_ATTCGAGG-GAGCTTGT_L003 032_AACCGTTC-TGAGCTAG_L003 \
003_TGAGCTAG-AACCGTTC_L004 033_ATATGCGC-CTGATCGT_L004 \
004_GAGACGAT-TAACCGGT_L003 034_CTTGGATG-GATAGCGA_L004 \
005_CTTGTCGA-GAACATCG_L003 035_TGTACCGT-TCTCTAGG_L004 \
006_TTCCAAGG-CCTTGTAG_L004 036_AAGTCGAG-ACATTGCG_L004 \
007_CGCATGAT-TCAGGCTT_L004 037_NNNNNNNN-NNNNNNNN_L003 \
008_ACGGAACA-GTTCTCGT_L004 038_CTACTTGG-CTGGAGTA_L003 \
009_ATCGATCG-TGCTTCCA_L003 039_NNNNNNNN-NNNNNNNN_L004 \
010_GCTATCCT-CACCTGTT_L003 040_ATCTCGCT-GCCTTGTT_L003 \
011_AGAGTAGC-TACGCCTT_L003 041_AGTTACGG-TTGGTCTC_L003 \
012_GACGATCT-ATGCACGA_L004 042_AGTTGGCT-ATGGTTGC_L004 \
013_AACTGAGC-CCTGATTG_L003 043_AATACGCG-CTATCGCA_L004 \
014_GTGCCATA-ACTAGGAG_L003 044_GTGGTGTT-CTTACAGC_L003 \
015_AAGCACTG-GTTGACCT_L004 045_NNNNNNNN-NNNNNNNN_L003 \
016_GTCGAAGA-CAATGTGG_L004 046_GTGAGCTT-TTGATCCG_L004 \
017_GATTACCG-CAGTCCAA_L004 047_NNNNNNNN-NNNNNNNN_L004 \
018_GTTGTTCG-GAAGTTGG_L003 048_AGGCATAG-TGCCTCTT_L004 \
019_CGGTTGTT-CATACCAC_L003 049_GCTGTAAG-GCTCTGTA_L004 \
020_AACTGGTG-TCACGTTC_L004 050_ATTGCGTG-GCCATAAC_L003 \
021_GAAGTTGG-GTTGTTCG_L003 051_TGATACGC-ACTCGTTG_L003 \
022_CGCAATCT-ATGCCTGT_L004 052_CGTCTTGT-TCTCGCAA_L004 \
023_GATGTGTG-CATGGCTA_L004 053_GATCCATG-TGGAGTTG_L004 \
024_GATTGCTC-GTGAAGTG_L004 054_GTTAAGGC-TCGAAGGT_L003 \
025_ACGTTCAG-GCACAACT_L003 055_ACGTCGTA-GTTCATGG_L004 \
026_CGTGTGTA-TTCGTTGG_L004 056_CTGCGTAT-CGAACTGT_L004 \
027_CCTGATTG-AACTGAGC_L004 057_ACGGTCTT-TCCGAGTT_L004 \
028_ATCACACG-TACGCTAC_L003 058_GATTGGAG-TTCTCTCG_L004 \
029_CACCTGTT-GCTATCCT_L004 059_TGTCCAGA-ATTCTGGC_L004 \
030_TCAGGCTT-CGCATGAT_L004 060_TTGACAGG-CAGTCTTC_L004

# submit the job for the defined samples to the queue
do
qsub $dir/ReCognitION_complete_analysis.sh "$sample"
done