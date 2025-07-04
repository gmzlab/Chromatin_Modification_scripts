#!/bin/bash

# PARAMETERS:
# Working directory.
WD=$1
# Genome folder absolute path.
GENOME=$2
# Number of processors.
THREADS=$3
# Sample name not including input or histone.
CHIP=$4
# Wether reads are paired end or not.
PAIRED_END=$5
# Maximum distance between peaks in regions.
MAX_DIST=$6
# Input name.
INPUT=$7
# Histone name.
HISTONE=$8

if [ $PAIRED_END == 'yes' ]
then

if [ $INPUT != 'not_provided' ] && [ $HISTONE != 'not_provided' ]
then

# Calling peaks on both input and histones but separately
epic2 --guess-bampe --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --control ${WD}/samples/${INPUT}/${INPUT}*.bam --chromsizes ${GENOME}/chromsizes > ${CHIP}_over_${INPUT}.peaks.bed

epic2 --guess-bampe --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --control ${WD}/samples/${HISTONE}/${HISTONE}*.bam --chromsizes ${GENOME}/chromsizes > ${CHIP}_over_${HISTONE}.peaks.bed

# Extracting high confidence sets by intersection of both peak calling modes (vs input and histone)
# -u: (unique) if any overlap exists, the A feature is reported.
# grep -v '_': excludes peaks on chromosomes called like chr1_..._random.
bedtools intersect -u -a ${CHIP}_over_${HISTONE}.peaks.bed -b ${CHIP}_over_${INPUT}.peaks.bed | grep -v '_' > ${CHIP}.peaks.bed

# Cleaning up.
rm ${CHIP}_over_*

elif [ $INPUT != 'not_provided' ]
then

# Calling peaks on input
epic2 --guess-bampe --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --control ${WD}/samples/${INPUT}/${INPUT}*.bam --chromsizes ${GENOME}/chromsizes | grep -v '_'  > ${CHIP}.peaks.bed

elif [ $HISTONE != 'not_provided' ]
then

# Calling peaks on histone
epic2 --guess-bampe --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --control ${WD}/samples/${HISTONE}/${HISTONE}*.bam --chromsizes ${GENOME}/chromsizes | grep -v '_'  > ${CHIP}.peaks.bed

else

# Calling peaks without control samples
epic2 --guess-bampe --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --chromsizes ${GENOME}/chromsizes | grep -v '_'  > ${CHIP}.peaks.bed

fi

else

if [ $INPUT != 'not_provided' ] && [ $HISTONE != 'not_provided' ]
then


# Extracting high confidence sets by intersection of both peak calling modes (vs input and histone)
# -u: (unique) if any overlap exists, the A feature is reported.
# grep -v '_': excludes peaks on chromosomes called like chr1_..._random.
bedtools intersect -u -a ${CHIP}_over_${HISTONE}.peaks.bed -b ${CHIP}_over_${INPUT}.peaks.bed | grep -v '_' > ${CHIP}.peaks.bed

# Cleaning up.
rm ${CHIP}_over_*.peaks.bed

elif [ $INPUT != 'not_provided' ]
then

# Calling peaks on input
epic2 --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --control ${WD}/samples/${INPUT}/${INPUT}*.bam --chromsizes ${GENOME}/chromsizes | grep -v '_' > ${CHIP}.peaks.bed

elif [ $HISTONE != 'not_provided' ]
then

# Calling peaks on histone
epic2 --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --control ${WD}/samples/${HISTONE}/${HISTONE}*.bam --chromsizes ${GENOME}/chromsizes | grep -v '_' > ${CHIP}.peaks.bed

else

# Calling peaks without control samples
epic2 --treatment ${WD}/samples/${CHIP}/${CHIP}*.bam --chromsizes ${GENOME}/chromsizes | grep -v '_' > ${CHIP}.peaks.bed

fi

fi

# Collapse contiguous peaks into single regions.
# -d: maximun distance between features allowed to be merged.
bedtools merge -d ${MAX_DIST} -i ${CHIP}.peaks.bed > ${CHIP}.regions.bed

# When analysing ChIP-seq data on human cells MAX_DIST was set to 5000 bp and when analysing data on yeast cells, it was set to 1000 bp.
