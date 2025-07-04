#!/bin/bash

if [ $# -ne 1 ]
then
	echo 'Parameters file not provided. Quiting...'
	exit
fi

############################# LOADING PARAMETERS #############################

echo ''
echo 'LOADING PARAMETERS'

# Working directory: absolute path to the folder containing all files.
WD=$(grep working_directory: $1 | awk '{ print $2 }')
# Genome folder absolute path.
GENOME=$(grep genome_path: $1 | awk '{ print $2 }')
# Number of processors.
THREADS=$(grep number_processors: $1 | awk '{ print $2 }')
# Wether reads are paired end or not.
PAIRED_END=$(grep paired_end: $1 | awk '{ print $2 }')
# Number of samples.
NSAMPLES=$(grep sample_ $1 | wc -l)
# Samples names saved within an array.
SAMPLES=()
for i in `seq 1 $NSAMPLES`
do
	SAMPLES[$i]=$(grep sample_${i}: $1 | awk '{ print $2 }')
done
# Number of replicates per sample, input and/or histone.
NREP=()
for i in `seq 1 $NSAMPLES`
do
	NREP[$i]=`ls --format single-column ${WD}/samples/${SAMPLES[$i]}/${SAMPLES[$i]}*.fastq.gz | wc -l`
	if [ $PAIRED_END == 'yes' ]
	then
		NREP[$i]=$((${NREP[$i]}/2))
	fi
done

############################# REPORTING PARAMETERS #############################

echo ''
echo ' - Working directory:' $WD
echo ' - Genome folder path:' $GENOME
echo ' - Number of processors:' $THREADS
echo ' - Paired end:' $PAIRED_END
echo ' - Number of samples:' $NSAMPLES
echo ' - Samples names:' ${SAMPLES[@]}
echo ' - Replicates per sample:' ${NREP[@]}
echo ''

sleep 10s

echo ''
echo 'PIPELINE RUNNING STARTED'
INIT=`date`
echo ''
echo ''

############################# BUILDING GENOME INDEX ############################# 

#cd $GENOME
#mkdir STAR_index
#STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles genome.fa


############################# READ MAPPING #############################

echo 'READ MAPPING'
echo ''

for i in `seq 1 $NSAMPLES`
do
	echo 'Processing' ${SAMPLES[$i]}':'
	echo ''
	cd ${WD}/samples/${SAMPLES[$i]}
	#fastqc ${SAMPLES[$i]}*.fastq*
	#mv ${SAMPLES[$i]}*fastqc* -t ../../quality_control
	#echo ''
	#echo ''
	echo 'Gunzipping...'
	gunzip -k ${SAMPLES[$i]}*.fastq.gz
	
	for j in `seq 1 ${NREP[$i]}`
	do
		if [ $PAIRED_END == 'yes' ]
		then
			STAR --runThreadN $THREADS --genomeDir $GENOME/STAR_index --readFilesIn ${SAMPLES[$i]}_${j}_1.fastq ${SAMPLES[$i]}_${j}_2.fastq --outSAMtype BAM Unsorted
		else
			STAR --runThreadN $THREADS --genomeDir $GENOME/STAR_index --readFilesIn ${SAMPLES[$i]}_${j}.fastq --outSAMtype BAM Unsorted
		fi
		# Filtering unique mapping reads
		# -F: exclude alignments with the given flag. FLAG 256 refers to secondary alignments.
		#samtools view -F 256 -O BAM -@ $THREADS Aligned.out.bam > ${SAMPLES[$i]}_smap.bam

		# USING MARKDUP.
		# Sorting .bam file by coordinates.
		# -@: number of threads.
		# -o: output file.
		samtools sort -@ $THREADS Aligned.out.bam -o ${SAMPLES[$i]}_${j}.bam 
	
		#mv Aligned.out.bam ${SAMPLES[$i]}.bam
	
		# Indexing .bam file
		samtools index -@ $THREADS ${SAMPLES[$i]}_${j}.bam
		# Generating bigwig files
		bamCoverage -p $THREADS -b ${SAMPLES[$i]}_${j}.bam --normalizeUsing RPKM -o ${SAMPLES[$i]}_${j}.bw 
		bamCoverage -p $THREADS -b ${SAMPLES[$i]}_${j}.bam --normalizeUsing RPKM -o ${SAMPLES[$i]}_${j}_fw.bw --filterRNAstrand forward
		bamCoverage -p $THREADS -b ${SAMPLES[$i]}_${j}.bam --normalizeUsing RPKM -o ${SAMPLES[$i]}_${j}_rv.bw --filterRNAstrand reverse
		#bigwigCompare -p $THREADS -b1 ${SAMPLES[$i]}_${j}_fw.bw -b2 ${SAMPLES[$i]}_${j}_rv.bw --operation add -o ${SAMPLES[$i]}_${j}_sum.bw
		
		# Cleaning up
		rm Aligned.out.bam ${SAMPLES[$i]}_${j}*.fastq
		
		if [ -d ${WD}/quality_control/${SAMPLES[$i]}_${j}_STAR_logs ]
		then
			echo 'Moving STAR logs to' ${WD}/quality_control/${SAMPLES[$i]}_${j}_STAR_logs
		else
			echo 'Creating' ${WD}/quality_control/${SAMPLES[$i]}_${j}_STAR_logs 'and moving STAR logs to it'
			mkdir ${WD}/quality_control/${SAMPLES[$i]}_${j}_STAR_logs
		fi
		mv Log.final.out  Log.out  Log.progress.out  SJ.out.tab -t ${WD}/quality_control/${SAMPLES[$i]}_${j}_STAR_logs
	done
	
	
	
	
	echo ''
	echo ''
	
done


# Relocating mapping files
cd ${WD}/samples
#mv */*.bam -t ../bam_files
#mv */*.bam.bai -t ../bam_files
#cp */*.bw -t ../bw_files
#mv */*.bw -t ../bw_files/per_replicate


echo 'PIPELINE RUNNING STARTED ON'
echo $INIT
echo ''
echo 'AND FINISHED ON'
date
echo ''

exit


###########################################################################

# WHEN USING ALL SAMPLES TOGETHER


	if [ $PAIRED_END == 'yes' ]
	then
	REPS1=${SAMPLES[$i]}_1_1.fastq
	REPS2=${SAMPLES[$i]}_1_2.fastq
	for j in `seq 2 ${NREP[$i]}`
	do
		REPS1=${REPS1},${SAMPLES[$i]}_${j}_1.fastq
		REPS2=${REPS2},${SAMPLES[$i]}_${j}_2.fastq
	done
	STAR --runThreadN $THREADS --genomeDir $GENOME/STAR_index --readFilesIn $REPS1 $REPS2 --outSAMtype BAM Unsorted
	
	else
	REPS=${SAMPLES[$i]}_1.fastq
	for j in `seq 2 ${NREP[$i]}`
	do
		REPS=${REPS},${SAMPLES[$i]}_${j}.fastq
	done
	STAR --runThreadN $THREADS --genomeDir $GENOME/STAR_index --readFilesIn $REPS --outSAMtype BAM Unsorted
	
	fi
	
	# Filtering unique mapping reads
	# -F: exclude alignments with the given flag. FLAG 256 refers to secondary alignments.
	#samtools view -F 256 -O BAM -@ $THREADS Aligned.out.bam > ${SAMPLES[$i]}_smap.bam

	# USING MARKDUP.
	# Sorting .bam file by coordinates.
	# -@: number of threads.
	# -o: output file.
	samtools sort -@ $THREADS Aligned.out.bam -o ${SAMPLES[$i]}.bam 
	
	#mv Aligned.out.bam ${SAMPLES[$i]}.bam
	
	# Indexing .bam file
	samtools index -@ $THREADS ${SAMPLES[$i]}.bam
	
	# Generating bigwig files
	bamCoverage -p ${THREADS} -b ${SAMPLES[$i]}.bam --normalizeUsing RPKM -o ${SAMPLES[$i]}.bw ###--filterRNAstrand
	
	# Cleaning up
	#rm Aligned.out.bam ${SAMPLES[$i]}*.fastq
	if [ -d ${WD}/quality_control/${SAMPLES[$i]}_STAR_logs ]
	then
		echo 'Moving STAR logs to' ${WD}/quality_control/${SAMPLES[$i]}_STAR_logs
	else
		echo 'Creating' ${WD}/quality_control/${SAMPLES[$i]}_STAR_logs 'and moving STAR logs to it'
		mkdir ${WD}/quality_control/${SAMPLES[$i]}_STAR_logs
	fi
	mv Log.final.out  Log.out  Log.progress.out  SJ.out.tab -t ${WD}/quality_control/${SAMPLES[$i]}_STAR_logs
	echo ''
	echo ''








