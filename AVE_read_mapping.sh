#!/bin/bash

# LOADING PARAMETERS (inherited from the parent function):

# Number of processors.
THREADS=$1
# Path to genome file.
GENOME_FILE=$2
# Sample name including input or histone.
SAMPLE=$3
# Wether reads are paired end or not.
PAIRED_END=$4
# Bin size for bigwig generation (number of bases).
BIN_SIZE=$5
# Wether there is a step for mapping to spike-in genome afterwards.
SPIKE=$6

#echo $THREADS $GENOME_FILE $SAMPLE $PAIRED_END

# MAPPING:

echo 'Gunzipping...'

if [ $PAIRED_END == 'yes' ]
then
	## PAIRED-END MAPPING
	
 	gunzip -k ${SAMPLE}_1.fastq.gz
	gunzip -k ${SAMPLE}_2.fastq.gz
	
	echo 'Number of reads in the original' $SAMPLE 'paired end .fastq files'
	expr `wc -l ${SAMPLE}_1.fastq | awk '{ print $1 }'` / 4
	expr `wc -l ${SAMPLE}_2.fastq | awk '{ print $1 }'` / 4
	echo ''
	

	# Mapping with BWA aligner:
	# -n: maximum edit distance.
	# -k: maximum edit distance in the seed.
	# -R: maximum number of equally best hits before proceeding with suboptimal alignments.
	# -t: number of processors.
	bwa aln -n 3 -k 2 -R 300 -t $THREADS $GENOME_FILE ${SAMPLE}_1.fastq > ${SAMPLE}_1.sai
	bwa aln -n 3 -k 2 -R 300 -t $THREADS $GENOME_FILE ${SAMPLE}_2.fastq > ${SAMPLE}_2.sai
	
 	# -n: maximum number of aligments to output in the XA (alternative hits) tag for reads paired properly.
	bwa sampe -n 3 $GENOME_FILE ${SAMPLE}_1.sai ${SAMPLE}_2.sai ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq > ${SAMPLE}.sam
	
	# Cleaning up:
	# -f: ignore if file does not exist.
	rm -f  ${SAMPLE}_1.sai ${SAMPLE}_2.sai 
	if [ $SPIKE != 'yes' ]
	then
		rm -f ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq
	fi
else
	
 	## SINGLE-END MAPPING
	
 	gunzip -k ${SAMPLE}.fastq.gz
	
	echo 'Number of reads in the original' $SAMPLE 'single end .fastq file'
	expr `wc -l ${SAMPLE}.fastq | awk '{ print $1 }'` / 4
	echo ''
	
	
	# Mapping with BWA aligner:
	# -n: maximum edit distance.
	# -k: maximum edit distance in the seed.
	# -R: maximum number of equally best hits before proceeding with suboptimal alignments.
	# -n: maximum number of aligments to output in the XA (alternative hits) tag for reads paired properly.
	bwa aln -n 3 -k 2 -R 300 -t $THREADS $GENOME_FILE ${SAMPLE}.fastq > ${SAMPLE}.sai
	bwa samse -n 3 $GENOME_FILE ${SAMPLE}.sai ${SAMPLE}.fastq > ${SAMPLE}.sam
	
	# Cleaning up:
	# -f: ignore if file does not exist.
	rm -f ${SAMPLE}.sai
	if [ $SPIKE != 'yes' ]
	then
		rm -f ${SAMPLE}.fastq
	fi

fi

# Filtering unique and high quality mapping reads:
# -F: exclude alignments with the given flag. FLAG 780 refers to unmapped reads, unmapped mates, secondary alignment and reads that fail platform/vendor quality checks (https://broadinstitute.github.io/picard/explain-flags.html).
# -q: exclude reads with mapping quality under a given value. 30 means a 0.001 probability of mismapping (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577856/).
# -O: output format.
# -@: number of processors.
samtools view -F 780 -O BAM -q 30 -@ $THREADS ${SAMPLE}.sam > ${SAMPLE}_smap.bam

# PCR DUPLICATES REMOVAL:

# Sorting .bam file by read name:
# -@: number of processors.
# -n: sort by read names.
# -o: output file.
samtools sort -@ $THREADS -n ${SAMPLE}_smap.bam -o ${SAMPLE}_sorted.bam 

# Filling in mate coordinates and inserting size fields:
# -@: number of processors.
# -m: add mate score tags used later by markdup.
samtools fixmate -@ $THREADS -m ${SAMPLE}_sorted.bam ${SAMPLE}_fixmate.bam

# Sorting again but by coordinates:
# -@: number of processors.
# -o: output file.
samtools sort -@ $THREADS ${SAMPLE}_fixmate.bam -o ${SAMPLE}_rmdup_sorted.bam 

# Removing PCR duplicates:
# -f: file on which save the statistics.
# -@: number of processors.
# -r: remove duplicate reads.
samtools markdup -f ${SAMPLE}_duplicate_remove_stats.txt -@ $THREADS -r ${SAMPLE}_rmdup_sorted.bam ${SAMPLE}.bam

mv ${SAMPLE}_duplicate_remove_stats.txt -t ../../quality_control

echo ''
echo 'Number of reads in original .sam file:'
samtools view -c ${SAMPLE}.sam 
echo 'Number of reads after PCR duplicates removal:'
samtools view -c ${SAMPLE}.bam
echo ''

# Indexing .bam file:
# -@: number of processors.
samtools index -@ $THREADS ${SAMPLE}.bam

# Generating bigwig files without Spike-in:
# -p: number of processors.
# -b: input BAM file.
# -bs: bin size for the bigwig (default 50).
# -normalizeUsing: normalizing method (RPKM, CPM, BPM...). CPM is set rather than RPKM as there is no gene length in ChIP-seq analysis. In this case, RPKM and CPM lead to the same results but in a different scale.
# -o: output BIGWIG file.
# -e: extend mate reads to consider the original DNA fragment. Only for paired end mode. 
if [ $PAIRED_END == 'yes' ]
then
	bamCoverage -p $THREADS -b ${SAMPLE}.bam -bs $BIN_SIZE --normalizeUsing CPM -o ${SAMPLE}.bw -e 
else
	bamCoverage -p $THREADS -b ${SAMPLE}.bam -bs $BIN_SIZE --normalizeUsing CPM -o ${SAMPLE}.bw
fi


# Cleaning up:
# -f: ignore if file does not exist.
rm -f ${SAMPLE}.sam ${SAMPLE}_smap.bam ${SAMPLE}_sorted.bam ${SAMPLE}_rmdup_sorted.bam ${SAMPLE}_fixmate.bam 


