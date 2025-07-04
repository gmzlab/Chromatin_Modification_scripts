#!/bin/bash


############################# LOADING PARAMETERS #############################

# Working directory: absolute path to the folder containing all files.
WD=$(grep working_directory: $1 | awk '{ print $2 }')
# Genome folder absolute path.
GENOME=$(grep genome_path: $1 | awk '{ print $2 }')
# Number of processors.
THREADS=$(grep number_processors: $1 | awk '{ print $2 }')



## A. Calculating coverage over the designated bins for later analysis ##

############################# GETTING GENOMIC BINS #############################
cd $GENOME

## -w binsize
bedtools makewindows -g hg38.genome -w 50000 > genome_50kb.bed

############################# COVERAGE OVER REGIONS #############################

multiBamSummary BED-file -p 12 --BED genome_50kb.bed -b ../bam_files/*.bam -o counts/50kb_rawCount.npz --outRawCounts counts/50kb_rawCount.txt

cat output_rawCount.txt | sort -k1,1 -k2,2n > output_rawCount.sort.txt

	
		

## B. Preparing and training the HMM model for chromatin states

############################# CONVERTING BAMS TO BEDS #############################

for i in `ls bam_files_merge/h3k*.bam bam_files_merge/input*.bam | cut -d '/' -f 2 | cut -d '.' -f 1` 
do
	echo $i
	bedtools bamtobed -i bam_files_merge/${i}.bam > bed_files/${i}.bed
done

############################# BINARIZING GENOME COVERAGE #############################
cd ${WD}

java -jar /home/gmzlab/opt/ChromHMM/ChromHMM.jar BinarizeBed -b 50000 /home/gmzlab/Documents/ricardo/genome/homo_sapiens_ucsc/chromsizes bed_files metadata/cell_mark_table_bed.tsv results/binarize

############################# CHROMATIN STATES #############################

# Provide number of states and reference genome directory
java -jar ../../opt/ChromHMM/ChromHMM.jar LearnModel -p 12 -b 50000 binarize_agustin_50kb/ model_50kb 4 hg38 

#############################  #############################



## C. Reordering the states or columns of a model without relearning the model, and outputs a model file and emission and transition tables with the states reordered ##

############################# COLUMN REORDER #############################

java -jar ../../opt/ChromHMM/ChromHMM.jar Reorder -f marks_order.txt model_agustin_10kb/model_4.txt model_10kb_reordered
java -jar ../../opt/ChromHMM/ChromHMM.jar Reorder -f marks_order.txt -o states_order.tsv -r model_t2t/RPE_5_segments.bed model_t2t_reordered_final/RPE_5_segments_reordered.bed model_t2t/model_5.txt model_t2t_reordered_final

#############################  #############################





## D. Obtain the correlation between the emission parameters of differente models and a reference ##

############################# CORRELATION BETWEEN MODELS #############################

java -jar ../../opt/ChromHMM/ChromHMM.jar CompareModels

#############################  #############################



## E. Take a learned model and binarized data, and output a segmentation. ##

############################# NEW SEGMENTATION #############################

java -jar ../../opt/ChromHMM/ChromHMM.jar MakeSegmentation -b 50000 model_5kb_reordered/model_5.txt binarize_agustin_50kb/ model_5kb_segmented_50kb

#############################  #############################


## F. Extract individual states from the segmentation file; and divide them into isolated bins for normalization ##

############################# STATES INDIVIDUAL ANALYSIS #############################

grep E1 model_5kb_segmented_50kb/RPE_5_50kb_segments.bed > 50kb_analysis/50kb_5_E1.bed ;
grep E2 model_5kb_segmented_50kb/RPE_5_50kb_segments.bed > 50kb_analysis/50kb_5_E2.bed ;
grep E3 model_5kb_segmented_50kb/RPE_5_50kb_segments.bed > 50kb_analysis/50kb_5_E3.bed ;
grep E4 model_5kb_segmented_50kb/RPE_5_50kb_segments.bed > 50kb_analysis/50kb_5_E4.bed ;
grep E5 model_5kb_segmented_50kb/RPE_5_50kb_segments.bed > 50kb_analysis/50kb_5_E5.bed

bedtools intersect -u -a ../../Data/genome/genome_50kb.bed -b 50kb_analysis/50kb_5_E1.bed > 50kb_analysis/E1_bins_50kb.bed ; 
bedtools intersect -u -a ../../Data/genome/genome_50kb.bed -b 50kb_analysis/50kb_5_E2.bed > 50kb_analysis/E2_bins_50kb.bed ;
bedtools intersect -u -a ../../Data/genome/genome_50kb.bed -b 50kb_analysis/50kb_5_E3.bed > 50kb_analysis/E3_bins_50kb.bed ;
bedtools intersect -u -a ../../Data/genome/genome_50kb.bed -b 50kb_analysis/50kb_5_E4.bed > 50kb_analysis/E4_bins_50kb.bed ;
bedtools intersect -u -a ../../Data/genome/genome_50kb.bed -b 50kb_analysis/50kb_5_E5.bed > 50kb_analysis/E5_bins_50kb.bed 

#############################  #############################




## G. Use the annotated genes instead of bins for intersect segmentation ##

########################## E1 GENES SCATTERPLOT SELECTION #######################
# -f: percentage of minimum overlap with target gene
bedtools intersect -u -f 0.5 -a Data/genome/Homo_sapiens_ann_genes_110_noChr.bed -b results/chromHMM/5kb_analysis/5kb_5_E1.bed > results/chromHMM/5kb_analysis/E1_genes.bed

multiBamSummary BED-file -p 12 --BED results/chromHMM/5kb_analysis/E1_genes.bed -b Data/mutK37R/h3k37me3_B12_1.bam Data/WT/bam_files_merge/h3k9me3.bam Data/mutK37R/input2_B12_1.bam Data/WT/bam_files_merge/input.bam -o results/chromHMM/5kb_analysis/genes_rawCount.npz --outRawCounts results/chromHMM/5kb_analysis/genes_rawCount.txt
