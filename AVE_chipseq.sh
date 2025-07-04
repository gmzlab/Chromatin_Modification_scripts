#!/bin/bash

if [ $# -lt 1 ]
then
	echo 'Parameters file not properly provided. Quiting...'
	exit
elif [ $# -eq 1 ]
then
	RUNMODE=$2dont
else
	RUNMODE=$2
fi


############################# LOADING PARAMETERS (I) #############################
# You must provide a .txt containing these parameters as the only argument when calling the AVE_chipseq.sh function.

# Working directory: absolute path to the folder containing all files.
WD=$(grep working_directory: $1 | awk '{ print $2 }')
# Genome folder absolute path.
GENOME=$(grep genome_path: $1 | awk '{ print $2 }')
# Wether Spike-in is included or not. Spiking a sample with DNA from a different organism and performing a secondary mapping allows normalization accounting for antibody effectiveness.
SPIKE=$(grep spikein: $1 | awk '{ print $2 }')
# Spike-in genome folder absolute path.
GENOME_SPIKE=$(grep genome_spike: $1 | awk '{ print $2 }')
# Number of processors.
THREADS=$(grep number_processors: $1 | awk '{ print $2 }')
# Wether reads are paired end or not.
PAIRED_END=$(grep paired_end: $1 | awk '{ print $2 }')
# Maximum distance between peaks in regions.
MAX_DIST=$(grep max_distance: $1 | awk '{ print $2 }')
# Bin size for bigwig generation (number of bases).
BIN_SIZE=$(grep bin_size: $1 | awk '{ print $2 }')
# Number of samples.
NSAMPLES=$(grep sample_ $1 | wc -l)
# Samples names saved within an array.
SAMPLES=()
for i in `seq 1 $NSAMPLES`
do
	SAMPLES[$i]=$(grep sample_${i}: $1 | awk '{ print $2 }')
done
# Input name.
INPUT=$(grep input: $1 | awk '{ print $2 }')
if [ -z $INPUT ]
then
	INPUT='not_provided'
fi
# Histone name.
HISTONE=$(grep histone: $1 | awk '{ print $2 }')
if [ -z $HISTONE ]
then
	HISTONE='not_provided'
fi



############################# WORKING DIRECTORY STRUCTURE ############################# 

if [ $RUNMODE == 'wd' ]
then
	echo 'Creating working directory on' $WD
	echo 'Please transfer the samples .fastq files to corresponding folders'
	echo 'and run this pipeline without setting the 2nd parameter to wd.'
	mkdir -p $WD
	cd $WD
	mkdir bam_files bam_files/per_replicate bw_files bw_files/per_replicate bw_files/spikein bw_files/normalized bw_files/normalized/spikein quality_control samples results results/plots results/clusters results/peaks
	cd samples
	mkdir ${SAMPLES[@]} 
	exit
fi


############################# LOADING PARAMETERS (II) #############################

# Number of replicates per sample, input and/or histone.
NREP=()
for i in `seq 1 $NSAMPLES`
do
	NREP[$i]=`ls --format single-column ${WD}/samples/${SAMPLES[$i]}/*.fastq.gz | wc -l`
	if [ $PAIRED_END == 'yes' ]
	then
		NREP[$i]=$((${NREP[$i]}/2))
	fi
done



############################# REPORTING PARAMETERS #############################

echo ''
echo 'LOADING PARAMETERS'

echo ''
echo ' - Working directory:' $WD
echo ' - Genome folder path:' $GENOME
echo ' - Spike-in normalization:' $SPIKE
echo ' - Spike-in genome folder path:' $GENOME_SPIKE
echo ' - Number of processors:' $THREADS
echo ' - Paired end:' $PAIRED_END
echo ' - Maximum distance between peaks:' $MAX_DIST
echo ' - Bin size for bigwig generation:' $BIN_SIZE
echo ' - Number of samples:' $NSAMPLES
echo ' - Samples names:' ${SAMPLES[@]}
echo ' - Replicates per sample:' ${NREP[@]}
echo ' - Input sample:' $INPUT
echo ' - Histone sample:' $HISTONE
if [ $RUNMODE != 'nopeak' ]
then
	echo ' - Run mode: read mapping and peak calling'
else
	echo ' - Run mode: read mapping only'
fi
echo ''

# Make sure that all the parameters are correct. If not, stop the pipeline and correct them.

sleep 10s


echo ''
echo 'PIPELINE RUNNING STARTED'
INIT=`date`
echo ''
echo ''



############################# BUILDING GENOME INDEX ############################# 

# This is a backup section for new assemblies when genome index is not built yet.

#cd $GENOME
#date
#bwa index *.fa
#samtools faidx *.fa
#cut -f1,2 *.fa.fai > chromsizes
#date

############################# READ MAPPING #############################

echo 'READ MAPPING'
echo ''
echo ''

for i in `seq 1 $NSAMPLES`
do
	cd ${WD}/samples/${SAMPLES[$i]}
	
	
	# QUALITY CONTROL:
	echo ${SAMPLES[$i]} 'quality control:'
	fastqc *.fastq.gz
	mv *fastqc* -t ../../quality_control
	echo ''
	echo ''
	
	
	# MAPPING EACH REPLICATE   --> calling the AVE_read_mapping.sh script
	for j in `seq 1 ${NREP[$i]}`
	do
		echo 'Processing' ${SAMPLES[$i]}_$j 'on'
		date
		echo ''
		read_mapping.sh $THREADS ${GENOME}/*.fa ${SAMPLES[$i]}_$j $PAIRED_END $BIN_SIZE $SPIKE
		
		if [ $SPIKE == 'yes' ]
		then
			read_mapping_spike.sh $THREADS ${GENOME_SPIKE}/*.fa ${SAMPLES[$i]}_${j} $PAIRED_END
			
			# When single end and spike the extend parameter will extend the reads to deepTools default parameter.
			bamCoverage -p $THREADS -b ${SAMPLES[$i]}_${j}.bam -o ${WD}/bw_files/spikein/${SAMPLES[$i]}_${j}.bw --binSize 500 -e 
			bamCoverage -p $THREADS -b ${SAMPLES[$i]}_${j}_spike.bam -o ${WD}/bw_files/spikein/${SAMPLES[$i]}_${j}_spike.bw --binSize 500 -e
			
			if [ ${SAMPLES[$i]} != $INPUT ]
			then
				cd ${WD}/bw_files/spikein
				# SPIKE-IN NORMALIZATION
				bigwigCompare -p $THREADS -b1 ${SAMPLES[$i]}_${j}.bw -b2 ${INPUT}_${j}.bw --skipZeroOverZero --skipNAs --operation ratio -o ${SAMPLES[$i]}_${j}_ratio.bedgraph -of bedgraph --fixedStep --binSize 500
				bigwigCompare -p $THREADS -b1 ${SAMPLES[$i]}_${j}_spike.bw -b2 ${INPUT}_${j}_spike.bw --skipZeroOverZero --skipNAs --operation ratio -o ${SAMPLES[$i]}_${j}_ratio_spike.bedgraph -of bedgraph --fixedStep --binSize 500
			
				ENRICH=$(echo "scale=4; `LC_NUMERIC="C" awk '{sum += $4} END {printf "%.5f\n", sum}' ${SAMPLES[$i]}_${j}_ratio.bedgraph`/`wc -l ${SAMPLES[$i]}_${j}_ratio.bedgraph | cut -f1 -d ' '`" | bc)
				ENRICH_SPIKE=$(echo "scale=4; `LC_NUMERIC="C" awk '{sum += $4} END {printf "%.5f\n", sum}' ${SAMPLES[$i]}_${j}_ratio_spike.bedgraph`/`wc -l ${SAMPLES[$i]}_${j}_ratio_spike.bedgraph | cut -f1 -d ' '`" | bc)
				echo ${SAMPLES[$i]}_${j}	`echo "scale=4; ${ENRICH}/${ENRICH_SPIKE}" | bc` >> scale_factors.txt
				
				cd ${WD}/samples/${SAMPLES[$i]}
				# When single end and spike the extend parameter will extend the reads to deepTools default parameter, and effectiveGenomeSize factor is added for human hg38.
				bamCoverage --bam ${SAMPLES[$i]}_${j}.bam -o ${WD}/bw_files/normalized/spikein/${SAMPLES[$i]}_${j}.bw --binSize $BIN_SIZE -e -p $THREADS --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --scaleFactor `grep $i ${WD}/bw_files/spikein/scale_factors.txt | cut -f2` 
			fi
		fi
		echo ''
	done
	
	
	# MERGING BAMS REPLICATES
	if [ ${NREP[$i]} -gt 1 ]
	then
		echo 'Merging' ${SAMPLES[$i]} 'bam files'
		samtools merge -@ $THREADS ${SAMPLES[$i]}_*.bam -o ${SAMPLES[$i]}.bam
		samtools index -@ $THREADS ${SAMPLES[$i]}.bam
		mv ${SAMPLES[$i]}.bam* -t ${WD}/bam_files
	else 
		echo 'There are no replicates to merge. Copying' ${SAMPLES[$i]} 'bam file to bam_files folder'
		cp ${SAMPLES[$i]}_1.bam ${WD}/bam_files/${SAMPLES[$i]}.bam
		cp ${SAMPLES[$i]}_1.bam.bai ${WD}/bam_files/${SAMPLES[$i]}.bam.bai
	fi
	
	
done



############################# PEAK CALLING #############################

# Calling AVE_peak_calling.sh if indicatedn so
if [ $RUNMODE != 'nopeak' ]
then
	echo 'PEAK CALLING'
	echo ''
	cd ${WD}/results

	for i in `seq 1 $NSAMPLES`
	do
	if [ ${SAMPLES[$i]} != $INPUT ] && [ ${SAMPLES[$i]} != $HISTONE ]
	then
		echo 'Processing' ${SAMPLES[$i]} 'on'
		date
		peak_calling.sh $WD $GENOME $THREADS ${SAMPLES[$i]} $PAIRED_END $MAX_DIST $INPUT $HISTONE 
		echo ''
		echo ''
	fi
	done

	mv *.peaks.* *.regions.* -t peaks
fi

# Relocating mapping files
echo 'Moving BAM files to' ${WD}/bam_files/per_replicate
cd ${WD}/samples
mv */*.bam -t ${WD}/bam_files/per_replicate
mv */*.bam.bai -t ${WD}/bam_files/per_replicate

#cd ${WD}/bam_files
#echo 'Sample,Total' > total_counts.csv
#for i in `seq 1 $NSAMPLES`
#do
#	for j in `seq 1 ${NREP[$i]}`
#	do
#		#echo ${SAMPLES[$i]}_$j,`samtools view -@ $THREADS -c ${SAMPLES[$i]}_$j.bam` >> total_counts.csv
#	done
#done

echo ''
echo ''



echo 'PIPELINE RUNNING STARTED ON'
echo $INIT
echo ''
echo 'AND FINISHED ON'
date
echo ''

exit

