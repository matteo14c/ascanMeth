#!/usr/bin/env bash

# Path to your modkit executable
MODKIT="/home/spini/dist_modkit_v0.5.0_5120ef7/modkit"

# Reference genome
REF_GENOME="/home/spini/ONTdata/hs1.fa"
REF_FAI="/home/spini/ONTdata/hs1.fa.fai"

# CPUs
CPUS=8


# Input BAM file
INPUT_BAM="MC_phased_SAMPLE_C.bam"
PREFIX="res"
M="${PREFIX}aScan_genome_M.csv"
C="${PREFIX}aScan_genome_C.csv"


# chromosomes

my_chr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22")

for chr in "${my_chr[@]}"; do
        OUT="${PREFIX}${chr}"
        echo Processing $chr
        $MODKIT pileup --region $chr --partition-tag HP --threads $CPUS --cpg --ref $REF_GENOME --combine-strands --filter-threshold C:0.75 --mod-threshold h:0.8 $INPUT_BAM $OUT 
        mv $OUT/1.bed $OUT/pat.bed
        mv $OUT/2.bed $OUT/mat.bed
	$MODKIT bm tobigwig --sizes $REF_FAI -m m $OUT/pat.bed $OUT/${chr}pat.bigWig
	$MODKIT bm tobigwig --sizes $REF_FAI -m m  $OUT/mat.bed $OUT/${chr}mat.bigWig
        cp aScanMeth.R $OUT
        cd $OUT
        Rscript aScanMeth.R --h1 pat.bed --h2 mat.bed --n 25 --c 10
	#rm aScanMeth.R
        cd ..
	cat $OUT/results_M_CpG.csv >> $M
	cat $OUT/results_D.CpG.csv >> $C
done
