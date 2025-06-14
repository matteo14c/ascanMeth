# Path to your modkit executable
MODKIT="../dist_modkit_v0.5.0_5120ef7/modkit"

# Reference genome
REF_GENOME="hs1.fa"



# Input BAM file
INPUT_BAM="phased_SAMPLE_P1.bam"
PREFIX="res_"

# chromosomes
my_chr=("chr4" "chr10")

for chr in "${my_chr[@]}"; do
	OUT="${PREFIX}${chr}"
	echo Processing $chr
	$MODKIT pileup --region $chr --partition-tag HP --cpg --ref $REF_GENOME --combine-strands --filter-threshold C:0.75 --mod-threshold h:0.8 $INPUT_BAM $OUT
	cp aScanMeth.R $OUT
        cd $OUT
	Rscript aScanMeth.R --h1 1.bed --h2 2.bed
	rm aScanMethMockUp.R
	cd ..	
done
