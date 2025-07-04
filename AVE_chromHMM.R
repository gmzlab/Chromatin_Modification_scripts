#################################### STATES - PIE CHART PLOT ############

library(GenomicRanges)

# Loading ChromHMM segmentation - 5kb / 5 states model
states.bed <- read.table('../results/chromHMM/model_agustin_5kb/RPE_5_segments.bed', sep = '\t') 
colnames(states.bed) <- c('chr', 'start', 'end', 'state')

# Removing extra scaffolds
states.bed <- states.bed[grepl("^.{1,2}$", states.bed$chr), ]
table(states.bed$chr)

# Creating GRanges object, and extracting length information for each segment
states <- GRanges(states.bed)

states$width <- width(GRanges(states.bed))-1
long <- c()
for (i in 1:5) {
  long[i] <- sum(states$width[states$state==paste0('E',i)])
}

long.100 <- round(long/sum(long)*100, 1)
long.100

pdf(file = 'plots/pie_chart_5kb.pdf', height = 4.5, width = 4.5)
par(lwd = 2)
pie(long, col = c('#9C9B9B','white','#2E519F', '#E19800', '#BB191A'), 
    main ='ChromHMM States', radius = 0.8, cex.main = 2 , border = "black",
    labels = c(paste0('S1 (',long.100[1],'%)'),
               paste0('S2 (',long.100[2],'%)'),
               paste0('S3 (',long.100[3],'%)'),
               paste0('S4 (',long.100[4],'%)'),
               paste0('S5 (',long.100[5],'%)')))
par(lwd = 1)
dev.off()


############################## LOADING AND PRE-PROCESSING COUNTS DATA ##############################

library(ggplot2)

# Loading raw counts table, for 14 marks + 3 inputs, quantified and compiled with deeptools multiBamSummary
tab <- read.csv2('../results/chromHMM/counts/5kb_RawCount.txt' , header = T, sep = '\t', quote ="'", check.names = T, as.is = T)
colnames(tab) <- gsub('.bam', '', colnames(tab))
colnames(tab)[1] <- 'chr'
rownames(tab) <- paste('bin', 1:nrow(tab), sep='')


# Removing extra scaffolds
tab <- tab[grepl("^.{1,2}$", tab$chr), ]

# Getting positions data.
bin.bed <- tab[,1:3]
bin.bed$binID <- rownames(bin.bed)

# Getting only count data as numeric.
marks.bin.cnt <- tab[,-(1:3)]
for (i in 1:ncol(marks.bin.cnt)) { marks.bin.cnt[,i] <- as.numeric(marks.bin.cnt[,i]) }
marks.bin.cnt <- as.matrix(marks.bin.cnt)



### DATA PROCESSING

# 1. CALCULATING CPM
marks.bin.cpm <- t(t(marks.bin.cnt)/(colSums(marks.bin.cnt)/1e6))

boxplot(marks.bin.cpm, ylim=c(0,7), main = "CPM", las = 2, cex.axis = 0.7, outline = FALSE)


# 2. REMOVING NULL ROWS
marks.bin.cpm.no0 <- marks.bin.cpm[rowSums(marks.bin.cpm) != 0, ]


# 3. Log2 + 1
marks.bin.log2 <- log2(marks.bin.cpm.no0 +1)

boxplot(marks.bin.log2, ylim=c(0,7), main = "log2", las = 2, cex.axis = 0.7, outline = FALSE)


# 4. INPUT CORRECTION
marks.corrected <- marks.bin.log2[,1:14]

# Reordering
marks.corrected <- marks.corrected[,c('h3k37me3_B12_1','h3k9me3','h3k27me3','h3k9me2','h3k9me1','h3k36me1','h3k36me2','h3k36me3','h3k37me1_B12','h3k4me1','h3k4me2','h3k4me3','h3k9ac','h3k27ac')]

marks.corrected[,c("h3k27ac", "h3k27me3","h3k36me1","h3k36me2","h3k36me3","h3k4me1",
  "h3k4me2","h3k4me3","h3k9ac","h3k9me1","h3k9me2","h3k9me3")] <-  marks.bin.log2[,c("h3k27ac", "h3k27me3","h3k36me1","h3k36me2","h3k36me3","h3k4me1",
                                                                                   "h3k4me2","h3k4me3","h3k9ac","h3k9me1","h3k9me2","h3k9me3")] - marks.bin.log2[, "input"]

marks.corrected[,"h3k37me1_B12"] <-  marks.bin.log2[,"h3k37me1_B12"] - marks.bin.log2[,"input1_B12"]
marks.corrected[,"h3k37me3_B12_1"] <-  marks.bin.log2[,"h3k37me3_B12_1"] - marks.bin.log2[,"input2_B12_1"]

# Renaming
colnames(marks.corrected) <-c('H3K37me3','H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2','H3K36me3','H3K37me1','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')

boxplot(marks.corrected, ylim=c(-1,7), main = "Input Corrected", las = 2, cex.axis = 0.7, outline = FALSE)


############################## LOADING EXTRA OMICS ##############################

omics <- read.csv2('../results/chromHMM/counts/5kb_extra_RawCount.txt' , header = T, sep = '\t', quote ="'", check.names = T, as.is = T)
colnames(omics)[1] <- 'chr'
rownames(omics) <- paste('bin', 1:nrow(omics), sep='')

omics <- omics[grepl("^.{1,2}$", omics$chr), ]


# Getting only count data as numeric.
omics.bin.cnt <- omics[,-(1:3)]
for (i in 1:ncol(omics.bin.cnt)) { omics.bin.cnt[,i] <- as.numeric(omics.bin.cnt[,i]) }
omics.bin.cnt <- as.matrix(omics.bin.cnt)

# Calculating CPM data.
omics.bin.cpm <- t(t(omics.bin.cnt)/(colSums(omics.bin.cnt)/1e6))

# Log2
omics.bin.log2 <- log2(omics.bin.cpm +1)


# INPUT CORRECTION

omics.corrected <- omics.bin.log2[,c("atac", "pol2", "h2az", "rna", "rna_nascent", "oris", "mcm3")]

omics.corrected[,"h2az"] <-  omics.bin.log2[,"h2az"] - omics.bin.log2[,"input2"]
omics.corrected[,"mcm3"] <-  omics.bin.log2[,"mcm3"] - omics.bin.log2[,"mcm3_input"]


boxplot(omics.corrected, ylim=c(-1,6), main = "Omics corrected", las = 2, cex.axis = 0.7, outline = F)


colnames(omics.corrected) <- c("ATAC", "Pol II", "H2A.Z", "RNA", "Nascent RNA", "Early Oris", "MCM3")

summary(omics.corrected)


############################## LOADING STATES BY 5Kb BINS ##############################

# These .bed had been generated by extracting the 5 states from RPE_5_segments.bed with grep, 
# and then bedtools intersect -u -a genome_5kb.bed -b ...5kb_5_E1.bed

E1.bed <- read.table('../results/chromHMM/model_agustin_5kb/E1_bins_5kb.bed', sep = '\t') 
colnames(E1.bed) <- c('chr', 'start', 'end')
E1.bed <- E1.bed[grepl("^.{1,2}$", E1.bed$chr), ]

E2.bed <- read.table('../results/chromHMM/model_agustin_5kb/E2_bins_5kb.bed', sep = '\t') 
colnames(E2.bed) <- c('chr', 'start', 'end')
E2.bed <- E2.bed[grepl("^.{1,2}$", E2.bed$chr), ]

E3.bed <- read.table('../results/chromHMM/model_agustin_5kb/E3_bins_5kb.bed', sep = '\t') 
colnames(E3.bed) <- c('chr', 'start', 'end')
E3.bed <- E3.bed[grepl("^.{1,2}$", E3.bed$chr), ]

E4.bed <- read.table('../results/chromHMM/model_agustin_5kb/E4_bins_5kb.bed', sep = '\t') 
colnames(E4.bed) <- c('chr', 'start', 'end')
E4.bed <- E4.bed[grepl("^.{1,2}$", E4.bed$chr), ]

E5.bed <- read.table('../results/chromHMM/model_agustin_5kb/E5_bins_5kb.bed', sep = '\t') 
colnames(E5.bed) <- c('chr', 'start', 'end')
E5.bed <- E5.bed[grepl("^.{1,2}$", E5.bed$chr), ]


# Turning data frames into data tables to make the merging easier
library(data.table)

setDT(states.bed)
setDT(E1.bed)
setDT(E2.bed)
setDT(E3.bed)
setDT(E4.bed)
setDT(E5.bed)

E1.bed[, state := "E1"]
E2.bed[, state := "E2"]
E3.bed[, state := "E3"]
E4.bed[, state := "E4"]
E5.bed[, state := "E5"]


# Merging all states tables, and merging this with the genomic bins table
all_E <- rbindlist(list(E1.bed, E2.bed, E3.bed, E4.bed, E5.bed), use.names = TRUE, fill = TRUE)

states.bed <- merge(bin.bed, all_E[, .(chr, start, end, state)],
                          by = c("chr", "start", "end"),
                          all.x = TRUE)

rownames(states.bed) <- states.bed$binID


# Extracting respective binIDs to states for later identification
E1.binID <- states.bed[states.bed$state %in% "E1", ]$binID
E2.binID <- states.bed[states.bed$state %in% "E2", ]$binID
E3.binID <- states.bed[states.bed$state %in% "E3", ]$binID
E4.binID <- states.bed[states.bed$state %in% "E4", ]$binID
E5.binID <- states.bed[states.bed$state %in% "E5", ]$binID

library(dplyr)
marks.enrich.state <- left_join(data.frame(marks.corrected, binID=rownames(marks.corrected)), states.bed, "binID")

omics.enrich.state <- left_join(data.frame(omics.corrected, binID=rownames(omics.corrected)), states.bed, "binID")

# Selecting corresponding rows in marks.bin.norm by matching rownames
marks.E1 <- marks.corrected[rownames(marks.corrected) %in% E1.binID,]
marks.E2 <- marks.corrected[rownames(marks.corrected) %in% E2.binID,]
marks.E3 <- marks.corrected[rownames(marks.corrected) %in% E3.binID,]
marks.E4 <- marks.corrected[rownames(marks.corrected) %in% E4.binID,]
marks.E5 <- marks.corrected[rownames(marks.corrected) %in% E5.binID,]

omics.E1 <- omics.corrected[rownames(omics.corrected) %in% E1.binID,]
omics.E2 <- omics.corrected[rownames(omics.corrected) %in% E2.binID,]
omics.E3 <- omics.corrected[rownames(omics.corrected) %in% E3.binID,]
omics.E4 <- omics.corrected[rownames(omics.corrected) %in% E4.binID,]
omics.E5 <- omics.corrected[rownames(omics.corrected) %in% E5.binID,]


## MERGING DATA INTO A MAIN MATRIX (with no E2 empty state and mitochondrial bins)

matrix <- left_join(marks.enrich.state, omics.enrich.state, 'binID')
matrix <- matrix[!is.na(matrix$state.x),]
matrix <- matrix[matrix$state.y != 'E2',]
matrix <- matrix[matrix$chr.x != 'MT',]

# Deleting all-null rows
matrix_no0 <- matrix[rowSums(matrix[,c(1:14,20:26)]) != 0, ]

save(matrix_no0, file = 'matrix_states.RData')


############################## GENERAL STATES BOXPLOTS ##############################
#load('matrix_states.RData')


# These graphs characterize the mark / omic distribution of each ChromHMM state

pdf(file = paste0('plots/boxplots_states_marks.pdf'), height = 15, width = 10)

colors <- c('#9C9B9B','#9C9B9B','#2E519F','#2E519F', '#BB191A', '#BB191A', '#E19800', '#BB191A','#BB191A','#BB191A','#BB191A','#BB191A','#BB191A','#BB191A')

par(mfrow=c(5,1))
boxplot(marks.E1, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,3), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "State 1", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:14, labels = colnames(marks.E1), tick =F, las = 2, cex.axis = 1.2, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

boxplot(marks.E2, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,3), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "State 2", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:14, labels = colnames(marks.E1), tick =F, las = 2, cex.axis = 1.2, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

boxplot(marks.E3, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,3), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "State 3", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:14, labels = colnames(marks.E1), tick =F, las = 2, cex.axis = 1.2, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

boxplot(marks.E4, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,3), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "State 4", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:14, labels = colnames(marks.E1), tick =F, las = 2, cex.axis = 1.2, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

boxplot(marks.E5, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,3), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "State 5", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:14, labels = colnames(marks.E1), tick =F, las = 2, cex.axis = 1.2, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

dev.off()







pdf(file = paste0('plots/boxplot_states_omics.pdf'), height = 15, width = 10)

colors <- colors <- c('#b543aa','#db5f1d','#db5f1d','#db5f1d', '#db5f1d', '#72b83d', '#72b83d')
par(mfrow=c(4,1))
boxplot(omics.E1, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,6), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "E1", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:7, labels = colnames(omics.E1), tick =F, las = 2, cex.axis = 0.9, hadj = 0.7)
abline(h=0, lwd=2, lty=2)


boxplot(omics.E3, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,6), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "E3", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:7, labels = colnames(omics.E1), tick =F, las = 2, cex.axis = 0.9, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

boxplot(omics.E4, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,6), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "E4", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:7, labels = colnames(omics.E1), tick =F, las = 2, cex.axis = 0.9, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

boxplot(omics.E5, boxlwd = 1.5 , lwd = 1.5 ,  ylim=c(-2.5,6), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
        col=colors, main = "E5", las = 2, cex.axis = 0.7, outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:7, labels = colnames(omics.E1), tick =F, las = 2, cex.axis = 0.9, hadj = 0.7)
abline(h=0, lwd=2, lty=2)

dev.off()


############################## INDIVIDUAL MARKS STATES ##############################
#load('matrix_states.RData')

# These graphs characterize each mark / omic throughout their ChromHMM states

## SELECTING ALL RELEVANT FEATURES
pdf(file = paste0('plots/boxplot_states_individual_all.pdf'), height = 6, width = 14)
colors <- c('#9C9B9B','#2E519F', '#E19800', '#BB191A')
par(mfrow=c(2,6))
marks <- c('H3K37me3','H3K9me3','H3K37me1','H3K36me3','H3K27me3','H3K36me2',"ATAC", "Pol.II","Nascent.RNA", "Oris")
mins <- c(-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,0,0,0,0)
maxs <- c(3,3,3,3,3,3,6,6,6,6)
for (i in 1:length(marks)) {
  mark1 <- marks[i]
  boxplot(matrix[[mark1]] ~ matrix$state.x, boxlwd = 1.5 , lwd = 1 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
          col=colors, main = mark1, las = 2, cex.axis = 0.7, outline = FALSE) 
  axis(2, cex.axis = 0.8)
  axis(1, at = 1:4, labels = c("CHt", "FHt", "Int", "Eu"), tick =F, las = 1, cex.axis = 1)
  abline(h=0, lwd=2, lty=2)
}
dev.off()



## ONLY SELECTED MARKS
pdf(file = paste0('plots/boxplot_states_individual_marks_clean.pdf'), height = 3, width = 12)
colors <- c('#9C9B9B','#2E519F', '#E19800', '#BB191A')
par(mfrow=c(1,5))
marks <- c('H3K37me3','H3K9me3','H3K27me3','H3K37me1', 'H3K36me3')
mins <- c(-2.5,-2.5,-2.5, -2.5, -2.5)
maxs <- c(3,3,3,3,3)
for (i in 1:length(marks)) {
  mark1 <- marks[i]
  boxplot(matrix[[mark1]] ~ matrix$state.x, boxlwd = 1.5 , lwd = 1 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "", 
          col=colors, main = mark1, las = 2, cex.main = 1.8, cex.axis = 0.7, outline = FALSE) 
  
  abline(h=0, lwd=2, lty=2)
}
dev.off()



## ONLY SELECTED OMICS
pdf(file = paste0('plots/boxplot_states_individual_omics_clean.pdf'), height = 3, width = 9)
colors <- c('#9C9B9B','#2E519F', '#E19800', '#BB191A')
par(mfrow=c(1,4))
marks <- c("ATAC", "Pol.II","Nascent.RNA", "Oris")
mins <- c(0,0,0,0)
maxs <- c(6,6,6,6)
for (i in 1:length(marks)) {
  mark1 <- marks[i]
  boxplot(matrix[[mark1]] ~ matrix$state.x, boxlwd = 1.5 , lwd = 1 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "log2(CPM+1)", 
          col=colors, main = gsub("\\.", " ", mark1), cex.main = 2 ,las = 2, cex.axis = 0.7, outline = FALSE) 
 
  abline(h=0, lwd=2, lty=2)
}
dev.off()





########################### SELECTED STATES ANALYSIS #######
#load('matrix_states.RData')


# Now we can select only one of the states for a dedicated analysis of the bins classifies as such.
matrix_no0.E1 <- matrix_no0[matrix_no0$state.x == 'E1',]

# We also select specific marks for comparison within this state

mark1 <-'H3K37me3'
mark2 <- 'H3K9me3'

matrix.selected <- matrix_no0.E1[,c(mark1, mark2, 'binID')] 

# Calculating bins differentially enriched for one of the marks
matrix.selected$ratio <- matrix.selected[,mark1]-matrix.selected[,mark2]
matrix.selected$class <- 'ns'
matrix.selected$class[matrix.selected$ratio < -1] <- 'k9'
matrix.selected$class[matrix.selected$ratio >  1] <- 'k37'

# Filtering matrixs
matrix_no0.E1 <-left_join(matrix_no0.E1,matrix.selected[, c('binID', 'class')],'binID' )
matrix_no0.E1.diff <- matrix_no0.E1[matrix_no0.E1$class != 'ns',]


###### FINAL SCATTERPLOT FOR ENRICHED BINS
volcol <- c( '#f77979', '#7171f5',"#9C9B9B")
names(volcol) <- c("k37","k9","ns")

pdf(paste0('plots/scatter_enrich_', mark1, '_', mark2, '.pdf'))
ggplot(matrix.selected, aes(x=H3K37me3 , y=H3K9me3, color=class)) +
  geom_point(size = 1) +
  scale_colour_manual(values = volcol) + 
  ylim(-2,3.1) + xlim(-2,3.1)+ 
  xlab(mark1) + ylab(mark2) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 19))+
  geom_abline(slope = 1, intercept = 1, color = "black", linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1, color = "black", linetype = "dashed") + 
  geom_abline(slope = lm(matrix.selected[, mark1] ~ matrix.selected[, mark2])$coefficients[2],
              intercept = lm(matrix.selected[, mark1] ~ matrix.selected[, mark2])$coefficients[1]
              , color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.7) + 
  annotate("text", x = -1.3, y = 2.7, label = paste0('cor = ', round(cor(matrix.selected[, mark1],matrix.selected[, mark2]),3)), size = 5, color = "red") +
  ggtitle(paste0(mark1, ' vs. ', mark2, ' within S1')) + 
  theme(plot.title = element_text(size = 18))

dev.off()



### BOXPLOTS FOR OMIC DATA IN ENRICHED BINS

# General
pdf(file = paste0('plots/boxplot_enriched_', mark1, '_', mark2, '.pdf'), height = 10, width = 10)
colors <- c( '#f77979', '#7171f5',"#9C9B9B")
par(mfrow=c(3,4))
omics.names <- c("ATAC", "Pol.II", "Nascent.RNA", "Oris",'H3K37me3','H3K9me3','H3K36me3', 'H3K37me1','H3K4me3','H3K27me3','H3K36me2','H3K9me1')
mins <- c(0,0,0,0, -2,-2,-1,-2, -1.5,-1.5,-3,-3)
maxs <- c(1.6,1.2,0.8,1.5, 3,3,1,1, 1, 1, 1, 1)
for (i in 1:length(omics.names)) {
  omic <- omics.names[i]
  boxplot(matrix_no0.E1.diff[[omic]] ~ matrix_no0.E1.diff$class, boxlwd = 2 , lwd = 2 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "log2(CPM)", 
          col=colors, main = omic, las = 2, cex.axis = 0.7, outline = FALSE) 
  axis(2, cex.axis = 0.8)
  axis(1, at = 1:2, labels =c('K37 enrich','K9 enrich'), tick =F, las = 1, cex.axis = 0.9)
  text(x=1.5, y =  maxs[i] * 1, labels = paste0("p = ", wilcox.test(matrix_no0.E1.diff[[omic]] ~ matrix_no0.E1.diff$class)$p.value), cex = 0.8)
  #abline(h=0, lwd=2, lty=2)
}

dev.off()


# Marks
pdf(file = paste0('plots/boxplot_enriched_marks_', mark1, '_', mark2, '.pdf'), height = 6, width = 7)
colors <- c( '#f77979', '#7171f5',"#9C9B9B")
par(mfrow=c(2,3))
omics.names <- c('H3K36me3','H3K4me3','H3K27me3','H3K36me2', 'H3K37me1','H3K9me1')
mins <- c(-1,-2, -1.5,-2,-2,-2)
maxs <- c(1,1, 1, 1, 1, 1)
for (i in 1:length(omics.names)) {
  omic <- omics.names[i]
  boxplot(matrix_no0.E1.diff[[omic]] ~ matrix_no0.E1.diff$class, boxlwd = 2 , lwd = 1.3 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "log2(CPM/Input)", 
          col=colors, main = omic, las = 2, cex.axis = 0.7, outline = FALSE) 
  axis(2, cex.axis = 0.8)
  axis(1, at = 1:2, labels =c('K37 enrich','K9 enrich'), tick =F, las = 1, cex.axis = 0.9)
  text(x=1.5, y =  maxs[i] * 1, labels = paste0("p = ", wilcox.test(matrix_no0.E1.diff[[omic]] ~ matrix_no0.E1.diff$class)$p.value), cex = 0.8)
  #abline(h=0, lwd=2, lty=2)
}
dev.off()


# General
pdf(file = paste0('plots/boxplot_enriched_omics_', mark1, '_', mark2, '.pdf'),height = 3, width = 15)
colors <- c( '#f77979', '#7171f5',"#9C9B9B")
par(mfrow=c(1,6))
omics.names <- c("ATAC", "Pol.II", "Nascent.RNA", "Oris")
mins <- c(0,0,0,0)
maxs <- c(1.6,1.2,0.8,1.5)
mains <- c('ATAC', 'Pol II', 'Nascent RNA', 'Early Oris')
for (i in 1:length(omics.names)) {
  omic <- omics.names[i]
  boxplot(matrix_no0.E1.diff[[omic]] ~ matrix_no0.E1.diff$class, boxlwd = 2 , lwd = 2 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "log2(CPM+1)", 
          col=colors, main = mains [i], las = 2, cex.axis = 0.7, outline = FALSE) 
  axis(2, cex.axis = 0.8)
  axis(1, at = 1:2, labels =c('K37 enrich','K9 enrich'), tick =F, las = 1, cex.axis = 0.9)
  text(x=1.5, y =  maxs[i] * 1, labels = paste0("p = ", wilcox.test(matrix_no0.E1.diff[[omic]] ~ matrix_no0.E1.diff$class)$p.value), cex = 0.8)
  #abline(h=0, lwd=2, lty=2)
}

dev.off()



################################ SELECTED STATES GENES & GO ########
#load('matrix_states.RData')

## Extracting enriched bins coordinates as a .bed file

matrix.enriched.k37me3.bed <- matrix_no0.E1.diff[matrix_no0.E1.diff$class == 'k37',c('chr.x', 'start.x', 'end.x')]
matrix.enriched.k9me3.bed <- matrix_no0.E1.diff[matrix_no0.E1.diff$class == 'k9',c('chr.x', 'start.x', 'end.x')]
matrix.enriched.ns.bed <- matrix_no0.E1[matrix_no0.E1$class == 'ns',c('chr.x', 'start.x', 'end.x')]


write.table(matrix.enriched.k37me3.bed, file = "../results/chromHMM/enriched_k37me3.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

 
write.table(matrix.enriched.k9me3.bed, file = "../results/chromHMM/enriched_k9me3.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(matrix.enriched.ns.bed, file = "../results/chromHMM/enriched_ns.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


write.table(matrix.enriched.ns.bed, file = "../results/chromHMM/enriched_ns.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


write.table(matrix_no0.E1[,c('chr.x', 'start.x', 'end.x','class')], file = "../results/chromHMM/5kb_analysis/E1_enriched.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Intersecting with genes.bed:
# bedtools intersect -u -f 0.1 -a ../../Data/genome/Homo_sapiens_ann_genes_110_noChr.bed -b enriched_k37me3.bed > genes_enriched_k37me3.bed

genes.enriched.k37  <- read.csv2('../results/chromHMM/5kb_analysis/genes_enriched_k37me3.bed' , header = F, sep = '\t', quote ="'", check.names = T, as.is = T)
genes.enriched.k9  <- read.csv2('../results/chromHMM/5kb_analysis/genes_enriched_k9me3.bed' , header = F, sep = '\t', quote ="'", check.names = T, as.is = T)

# Refined gene selection must cover at least 10% of the target gene
genes.enriched.k37.refined <- read.csv2('../results/chromHMM/5kb_analysis/genes_enriched_k37me3_refined.bed' , header = F, sep = '\t', quote ="'", check.names = T, as.is = T)
genes.enriched.k9.refined <- read.csv2('../results/chromHMM/5kb_analysis/genes_enriched_k9me3_refined.bed' , header = F, sep = '\t', quote ="'", check.names = T, as.is = T)

genes.enriched.k37.ID <- genes.enriched.k37[,4]
genes.enriched.k9.ID <- genes.enriched.k9[,4]

genes.enriched.k37.refined.ID <- genes.enriched.k37.refined[,4]
genes.enriched.k9.refined.ID <- genes.enriched.k9.refined[,4]


## GENE ONTOLOGY ANALYSIS
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(enrichplot))


###
go.enrich.k37 <- enrichGO(gene = genes.enriched.k37.ID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = 'BP',
                      pAdjustMethod = 'BH',
                      pvalueCutoff  = 0.05,
                      keyType = 'ENSEMBL')

barplot(go.enrich.k37, showCategory = 10, title = 'K37me3 enriched genes')
dotplot(go.enrich.k37, showCategory = 15, title = 'K37me3 enriched genes')

pdf(file = 'plots/go_enrich_k37.pdf', height = 7, width = 9)
dotplot(go.enrich.k37, showCategory = 15, title = 'K37me3 enriched genes')
dev.off()


###

go.enrich.k37.refined <- enrichGO(gene = genes.enriched.k37.refined.ID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = 'BP',
                          pAdjustMethod = 'BH',
                          qvalueCutoff  = 0.05,
                          keyType = 'ENSEMBL')

barplot(go.enrich.k37, showCategory = 10, title = 'K37me3 enriched genes')
dotplot(go.enrich.k37.refined, showCategory = 15, title = 'K37me3 enriched refined genes')

pdf(file = 'plots/go_enrich_k37_refined.pdf', height = 5, width = 6)
dotplot(go.enrich.k37.refined, showCategory = 15, title = 'K37me3 enriched genes')
dev.off()




###

go.enrich.k9 <- enrichGO(gene = genes.k9.test,
                          OrgDb         = org.Hs.eg.db,
                          ont           = 'BP',
                          pAdjustMethod = 'BH',
                          pvalueCutoff  = 0.05,
                          keyType = 'ENSEMBL')


barplot(go.enrich.k9, showCategory = 10, title = 'K9me3 enriched genes')
dotplot(go.enrich.k9, showCategory = 15, title = 'K9me3 enriched genes')



#### 


go.enrich.k9.refined <- enrichGO(gene = genes.enriched.k9.refined.ID,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = 'BP',
                                  pAdjustMethod = 'BH',
                                  pvalueCutoff  = 1,
                                  keyType = 'ENSEMBL')

barplot(go.enrich.k9.refined, showCategory = 10, title = 'K37me3 enriched genes')
dotplot(go.enrich.k9.refined, showCategory = 15, title = 'K37me3 enriched genes')

pdf(file = 'plots/go_enrich_k37.pdf', height = 7, width = 9)
dotplot(go.enrich.k9.refined, showCategory = 15, title = 'K37me3 enriched genes')


## COMMON GENES INTERSECTION

length(genes.enriched.k37.ID)
length(genes.enriched.k9.ID)
length(intersect(genes.enriched.k37.ID, genes.enriched.k9.ID))


length(genes.enriched.k37.refined.ID)
length(genes.enriched.k9.refined.ID)
length(intersect(genes.enriched.k37.refined.ID, genes.enriched.k9.refined.ID))


genes.k9.test <- genes.enriched.k9.ID[!genes.enriched.k9.ID %in% genes.enriched.k9.refined.ID]


library(VennDiagram)



pdf("plots/venn_genes.pdf", width = 2, height = 2)
venn.plot <- venn.diagram(
  x = list(H3K37me3 = genes.enriched.k37.refined.ID, H3K9me3 = genes.enriched.k9.refined.ID),
  filename = NULL, main = 'Enriched genes',            
  fill = c("red","blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main.cex = 1.6,
  cat.pos = 0)
grid.draw(venn.plot)                           
dev.off() 

### GO TERMS INTERSECTION
# Extraer ID de términos GO
terms.k37 <- go.enrich.k37$ID
terms.k9  <- go.enrich.k9$ID

length(go.enrich.k37$ID)
length(go.enrich.k9$ID)
length(intersect(go.enrich.k37$ID, go.enrich.k9$ID))


pdf("plots/venn_terms.pdf", width = 2, height = 2)
venn.plot <- venn.diagram(
  x = list(H3K37me3 = go.enrich.k37$ID, H3K9me3 = go.enrich.k9$ID),
  filename = NULL, main = 'Enriched GO terms',            
  fill = c("red","blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main.cex = 1.6,
  cat.pos = 0)
grid.draw(venn.plot)                           
dev.off() 



venn.plot <- venn.diagram(
  x = list(H3K37me3 = binarized.E1$binID[binarized.E1$H3K37me3 == 1], H3K9me3 = binarized.E1$binID[binarized.E1$H3K9me3 == 1]),
  filename = NULL, main = 'Enriched genes',            
  fill = c("red","blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main.cex = 1.6,
  cat.pos = 0)
grid.draw(venn.plot)   



####################################### E1 GENES SIGNAL ANALYSIS ########
# We now compare if the differentially enriched genes actually have a different signal for the marks for validation purposes

tab.genes <- read.csv2('../results/chromHMM/5kb_analysis/genes_rawCount_sorted.txt' , header = T, sep = '\t', quote ="'", check.names = T, as.is = T)
tab.genes.metadata <- read.csv2('../results/chromHMM/5kb_analysis/E1_genes.bed' , header = F, sep = '\t', quote ="'", check.names = T, as.is = T)
tab.genes$ID <- tab.genes.metadata[,4]

colnames(tab.genes) <- gsub('.bam', '', colnames(tab.genes))
colnames(tab.genes)[1] <- 'chr'
rownames(tab.genes) <- paste('bin', 1:nrow(tab.genes), sep='')





tab.genes.quant <- tab.genes[,4:7]
for (i in 1:ncol(tab.genes.quant)) { tab.genes.quant[,i] <- as.numeric(tab.genes.quant[,i]) }
tab.genes.quant <- as.matrix(tab.genes.quant)


boxplot(tab.genes.quant, las = 2, cex.axis = 0.7, outline = FALSE)

# Calculating CPM data.
tab.genes.quant <- t(t(tab.genes.quant)/(colSums(tab.genes.quant)/1e6))


#tab.genes.quant <- tab.genes.quant[rowSums(tab.genes.quant) != 0, ]


tab.genes.quant <- log2(tab.genes.quant +1)



# INPUT CORRECTION
tab.genes.quant.corrected <- tab.genes.quant[,1:2]


tab.genes.quant.corrected[,"h3k37me3_B12_1"] <-  tab.genes.quant[,"h3k37me3_B12_1"] - tab.genes.quant[,"input2_B12_1"]
tab.genes.quant.corrected[,"h3k9me3"] <-  tab.genes.quant[,"h3k9me3"] - tab.genes.quant[,"input"]



boxplot(tab.genes.quant.corrected, main = "Input Corrected", las = 2, cex.axis = 0.7, outline = FALSE)

# Selected marks

mark1 <-'h3k37me3_B12_1'
mark2 <- 'h3k9me3'

# Matrix of selected marks within the selected state (E1) 
matrix.selected <- as.data.frame(tab.genes.quant.corrected[,1:2])

# Calculating bins differentially enriched for one of the marks
matrix.selected$ratio <- matrix.selected[,mark1]-matrix.selected[,mark2]
matrix.selected$class <- 'ns'
matrix.selected$class[matrix.selected$ratio < -1] <- 'k9'
matrix.selected$class[matrix.selected$ratio >  1] <- 'k37'
matrix.selected$ID <- tab.genes$ID
matrix.selected <- cbind(matrix.selected, tab.genes [,1:3],tab.genes$ID)
table(matrix.selected$class)


###### FINAL SCATTERPLOT FOR ENRICHED BINS
volcol <- c( '#f77979', '#7171f5',"#9C9B9B")
names(volcol) <- c("k37","k9","ns")

pdf(paste0('plots/scatter_genes_', mark1, '_', mark2, '.pdf'))
ggplot(matrix.selected, aes(x=h3k37me3_B12_1 , y=h3k9me3, color=class)) +
  geom_point(size = 1) +
  scale_colour_manual(values = volcol) + 
 # ylim(-2,3.1) + xlim(-2,3.1)+ 
  xlab(mark1) + ylab(mark2) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 19))+
  geom_abline(slope = 1, intercept = 1, color = "black", linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1, color = "black", linetype = "dashed") + 
  geom_abline(slope = lm(matrix.selected[, mark1] ~ matrix.selected[, mark2])$coefficients[2],
              intercept = lm(matrix.selected[, mark1] ~ matrix.selected[, mark2])$coefficients[1]
              , color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.7) + 
  annotate("text", x = -1.3, y = 2.7, label = paste0('cor = ', round(cor(matrix.selected[, mark1],matrix.selected[, mark2]),3)), size = 5, color = "red") +
  ggtitle(paste0(mark1, ' vs. ', mark2, ' within S1')) + 
  theme(plot.title = element_text(size = 18))

dev.off()


##### GENES


genes.full.k37.ID <- matrix.selected$ID[matrix.selected$class == "k37"]
genes.full.k9.ID <- matrix.selected$ID[matrix.selected$class == "k9"]

length(genes.full.k37.ID)
length(genes.full.k9.ID)

length(intersect(genes.full.k37.ID, genes.full.k9.ID))

length(intersect(genes.full.k37.ID, genes.enriched.k37.hg38.ID))
length(intersect(genes.full.k37.ID, genes.enriched.k37.t2t.ENS))

length(intersect(genes.full.k9.ID, genes.enriched.k9.hg38.ID))
length(intersect(genes.full.k9.ID, genes.enriched.k9.t2t.ENS))




go.enrich.k37 <- enrichGO(gene = genes.full.k37.ID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = 'BP',
                          pAdjustMethod = 'BH',
                          pvalueCutoff  = 0.3,
                          keyType = 'ENSEMBL')





go.enrich.k9 <- enrichGO(gene = genes.full.k9.ID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = 'BP',
                          pAdjustMethod = 'BH',
                          pvalueCutoff  = 0.05,
                          keyType = 'ENSEMBL')


pdf(file = 'plots/go_enrich_k37.pdf', height = 7, width = 9)
dotplot(go.enrich.k37, showCategory = 15, title = 'K37me3 enriched genes')
dev.off()



### Saving BED for metageneplotting

matrix.genes.k37me3.bed <- matrix.selected[matrix.selected$class == 'k37',c('X.chr', 'start', 'end','tab.genes$ID')]
matrix.genes.k9me3.bed <- matrix.selected[matrix.selected$class == 'k9',c('X.chr', 'start', 'end', 'tab.genes$ID')]
matrix.genes.ns.bed <- matrix.selected[matrix.selected$class == 'ns',c('X.chr', 'start', 'end', 'tab.genes$ID')]


write.table(matrix.genes.k37me3.bed, file = "../results/metaplot/genes_k37me3.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(matrix.genes.k9me3.bed, file = "../results/metaplot/genes_k9me3.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(matrix.genes.ns.bed, file = "../results/metaplot/genes_ns.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





############### BINARIZED MARKS ANALYSIS #############

# Extracting binarized signals of selected marks in the cannonical chromosomes

extract_binarized_chromosomes <- function(path, columns = c("H3K37me3", "H3K9me3")) {
  files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)
  
  all_data <- lapply(files, function(file) {
    raw <- read.table(file, sep = "\t", header = FALSE, skip = 1, stringsAsFactors = FALSE)
    colnames(raw) <- raw[1, ]
    df <- raw[-1, columns, drop = FALSE]
    
    # Making sure that all values are numeric
    df[] <- lapply(df, as.numeric)
    
    # Extrating file ID
    filename <- basename(file)
    chrom <- sub("^.*RPE_([0-9XYM]+)_binary\\.txt$", "\\1", filename)
    df$chr <- chrom
    
    # Adding an additional row of 0s to account for the extra final bin
    zero_row <- data.frame(matrix(0, ncol = length(columns), nrow = 1))
    colnames(zero_row) <- columns
    zero_row$chr <- chrom
    
    df <- rbind(df, zero_row)
    return(df)
  })
  
  combined <- do.call(rbind, all_data)
  return(combined)
}


binarized.marks <- extract_binarized_chromosomes(path = '../results/chromHMM/binarize_agustin_5kb/')


table(binarized.marks$chr)

bin.bed <- bin.bed[bin.bed$chr != 'MT',]
binarized.marks <- cbind(binarized.marks, bin.bed[,2:4])

binarized.E1 <- binarized.marks[binarized.marks$binID %in% E1.binID, ]
#save(binarized.E1, file = 'binarized_E1.RData')

load('binarized_E1.RData')


# Extracting beds according to marks presence/absence

sum(binarized.E1$H3K37me3 == 1 & binarized.E1$H3K9me3 == 0)
sum(binarized.E1$H3K37me3 == 0 & binarized.E1$H3K9me3 == 1)
sum(binarized.E1$H3K37me3 == 1 & binarized.E1$H3K9me3 == 1)
sum(binarized.E1$H3K37me3 == 0 & binarized.E1$H3K9me3 == 0)


binarized.E1.k37me3 <- binarized.E1[binarized.E1$H3K37me3 == 1 & binarized.E1$H3K9me3 == 0, 3:6]
binarized.E1.k9me3 <- binarized.E1[binarized.E1$H3K37me3 == 0 & binarized.E1$H3K9me3 == 1, 3:6]
binarized.E1.both <- binarized.E1[binarized.E1$H3K37me3 == 1 & binarized.E1$H3K9me3 == 1, 3:6]
binarized.E1.none <- binarized.E1[binarized.E1$H3K37me3 == 0 & binarized.E1$H3K9me3 == 0, 3:6]

write.table(binarized.E1.k37me3, file = "../results/chromHMM/binarized_E1/binarized_k37me3test.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(binarized.E1.k9me3, file = "../results/chromHMM/binarized_E1/binarized_k9me3.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(binarized.E1.both, file = "../results/chromHMM/binarized_E1/binarized_both.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(binarized.E1.none, file = "../results/chromHMM/binarized_E1/binarized_none.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Combining binarized states with the main data matrix
matrix.binarized.E1 <- matrix_no0.E1
matrix.binarized.E1$binarized <- 'null'

matrix.binarized.E1$binarized[matrix.binarized.E1$binID %in% binarized.E1.k37me3$binID] <- "k37me3"
matrix.binarized.E1$binarized[matrix.binarized.E1$binID %in% binarized.E1.k9me3$binID] <- "k9me3"
matrix.binarized.E1$binarized[matrix.binarized.E1$binID %in% binarized.E1.both$binID] <- "both"
matrix.binarized.E1$binarized[matrix.binarized.E1$binID %in% binarized.E1.none$binID] <- "none"

table(matrix.binarized.E1$binarized)





###### SCATTERPLOT FOR BINARIZED BINS
volcol <- c( '#f77979', '#7171f5','#30bf56','#2d2e2d')
names(volcol) <- c("k37me3","k9me3","both", "none")

pdf(paste0('plots_binarized/scatter_binarized_', mark1, '_', mark2, '.pdf'))
ggplot(matrix.binarized.E1, aes(x=H3K37me3 , y=H3K9me3, color=binarized)) +
  geom_point(size = 1) +
  scale_colour_manual(values = volcol) + 
  ylim(-2,3.1) + xlim(-2,3.1)+ 
  xlab(mark1) + ylab(mark2) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 19))+
  geom_abline(slope = 1, intercept = 1, color = "black", linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1, color = "black", linetype = "dashed") + 
  geom_abline(slope = 0, intercept = 0, color = "black", linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.3, alpha = 0.7) +
  annotate("text", x = -1.3, y = 2.7, label = paste0('cor = ', round(cor(matrix.selected[, mark1],matrix.selected[, mark2]),3)), size = 5, color = "red") +
  ggtitle(paste0(mark1, ' vs. ', mark2, ' binarized in E1')) + 
  theme(plot.title = element_text(size = 18))

dev.off()



## BOXPLOTS FOR OMIC DATA IN ENRICHED BINS


levels(as.factor(matrix.binarized.E1$binarized))

matrix.binarized.E1$binarized <- factor(matrix.binarized.E1$binarized,levels = c("k37me3", "k9me3", "both", "none"))


pdf(file = paste0('plots_binarized/boxplot_binarized_', mark1, '_', mark2, '.pdf'), height = 10, width = 12)
colors <- c( '#f77979', '#7171f5','#30bf56','#2d2e2d')
par(mfrow=c(3,4))
feature.names <- c("ATAC", "Pol.II", "Nascent.RNA", "Oris",'H3K37me3','H3K9me3','H3K36me3', 'H3K37me1','H3K4me3','H3K27me3','H3K36me2','H3K9me1')
mins <- c(0,0,0,0, -2,-2,-1,-2, -1.5,-1.5,-3,-3)
maxs <- c(1.6,1.2,0.8,1.5, 3,3,1,1, 1, 1, 1, 1)

for (i in 1:length(feature.names)) {
  feature <- feature.names[i]
  boxplot(matrix.binarized.E1[[feature]] ~ matrix.binarized.E1$binarized, boxlwd = 2 , lwd = 2 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "log2(CPM)", 
          col=colors, main = feature, las = 2, cex.axis = 0.7, outline = FALSE) 
  axis(2, cex.axis = 0.8)
  axis(1, at = 1:4, labels =c('K37 only','K9 only', 'Both', 'None'), tick =F, las = 1, cex.axis = 0.9)
  #text(x=1.5, y =  maxs[i] * 1, labels = paste0("p = ", wilcox.test(matrix.binarized.E1[[feature]] ~ matrix.binarized.E1$binarized)$p.value), cex = 0.8)
  #abline(h=0, lwd=2, lty=2)
}
dev.off()




binarized.E1$binID[binarized.E1$H3K37me3 == 1]

sum(binarized.E1$H3K37me3 == 1 & binarized.E1$H3K9me3 == 0)
sum(binarized.E1$H3K37me3 == 0 & binarized.E1$H3K9me3 == 1)





############################### SMOOTHED E1 ANALYSIS #################

# Same analysis as before, but smoothing across a 50 kb window to reduce potential local noise

# Smoothing .bam files with: 
# bamCoverage -p 12 -b FILE.bam -bs 50 --normalizeUsing CPM --smoothLength 50000 -o FILE.bw"

# Quantifying .bw files in E1 5kb bins with: 
# multiBigwigSummary BED-file --BED .../E1_bins_5kb.bed -b *.bw -o multiBigWig_smoothed_50kb.npz --outRawCounts multiBigWig_smoothed_50kb.tsv -p 12 --smartLabels


smooth.E1.tsv <- read_tsv("../Data/WT/bw_files_merge/bw_smoothed_5kb_10kb/multiBigWig_E1_smoothed_5kb_10kb.tsv")

# Loading and merging omics
smooth.E1.omics.tsv <- read_tsv("../Data/WT/bw_files_merge/bw_smoothed_50000/multiBigWig_E1_omics_smoothed_50kb.tsv")
smooth.E1.tsv <- cbind(smooth.E1.tsv, smooth.E1.omics.tsv[,4:7])

# Cleaning
colnames(smooth.E1.tsv) <- gsub("'", "", colnames(smooth.E1.tsv))
smooth.E1.tsv <- smooth.E1.tsv[grepl("^.{1,2}$", smooth.E1.tsv$`#chr`), ]
smooth.E1.tsv <- smooth.E1.tsv[smooth.E1.tsv$`#chr` != 'MT', ]

# ??
#table(smooth.E1.tsv$`#chr`)
#table(matrix_no0.E1$chr.x)

# Adding (arbitrary) identifiers at the 4th column
smooth.E1.tsv <- cbind(
  smooth.E1.tsv[, 1:3],
  binID =  paste0("bin", seq_len(nrow(smooth.E1.tsv))),
  smooth.E1.tsv[, 4:ncol(smooth.E1.tsv)])

rownames(smooth.E1.tsv) <- smooth.E1.tsv$binID
smooth.bed <- smooth.E1.tsv[,1:4]
smooth.E1 <- smooth.E1.tsv[,-(1:4)]

### TRANSFORMING (from CPM) ###

smooth.E1 <- subset(smooth.E1, select = -h3k37me3_B10_1)
smooth.E1 <- as.matrix(smooth.E1)

# boxplot(smooth.E1, ylim=c(0,0.3), main = "CPM", las = 2, cex.axis = 0.7, outline = FALSE)

# No0
smooth.E1 <- smooth.E1[rowSums(smooth.E1) != 0, ]

## Log2 + 1
smooth.E1 <- log2(smooth.E1 +1)

# Input Correction
smooth.E1.corrected <- smooth.E1[,c(1:14,18:21)]
smooth.E1.corrected <- smooth.E1.corrected[,c('h3k37me3_B12_1','h3k9me3','h3k27me3','h3k9me2','h3k9me1','h3k36me1','h3k36me2','h3k36me3','h3k37me1_B12','h3k4me1','h3k4me2','h3k4me3','h3k9ac','h3k27ac','atac','oris','pol2','rna_nascent')]

smooth.E1.corrected[,c("h3k27ac", "h3k27me3","h3k36me1","h3k36me2","h3k36me3","h3k4me1",
                   "h3k4me2","h3k4me3","h3k9ac","h3k9me1","h3k9me2","h3k9me3")] <-  smooth.E1[,c("h3k27ac", "h3k27me3","h3k36me1","h3k36me2","h3k36me3","h3k4me1",
                                                                                                      "h3k4me2","h3k4me3","h3k9ac","h3k9me1","h3k9me2","h3k9me3")] - smooth.E1[, "input"]
smooth.E1.corrected[,"h3k37me1_B12"] <-  smooth.E1[,"h3k37me1_B12"] - smooth.E1[,"input1_B12"]
smooth.E1.corrected[,"h3k37me3_B12_1"] <-  smooth.E1[,"h3k37me3_B12_1"] - smooth.E1[,"input2_B12_1"]

colnames(smooth.E1.corrected) <-c('H3K37me3','H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2','H3K36me3','H3K37me1','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac','ATAC','Oris','PolII','Nascent RNA')

# boxplot(smooth.E1.corrected, ylim=c(-2,3), main = "Input Corrected", las = 2, cex.axis = 0.7, outline = FALSE)


smooth.E1.matrix <- cbind(smooth.bed, smooth.E1.corrected)


#### ANALYSIS

mark1 <-'H3K37me3'
mark2 <- 'H3K9me3'

# Matrix of selected marks within the selected state (E1) 
matrix.smooth.selected <- smooth.E1.matrix[,c(mark1, mark2, 'binID')] 

# Calculating bins differentially enriched for one of the marks
# NOW WE TAKE QUANTILES
matrix.smooth.selected$ratio <- matrix.smooth.selected[,mark1]-matrix.smooth.selected[,mark2]
matrix.smooth.selected$class <- 'ns'
matrix.smooth.selected$class[matrix.smooth.selected$ratio < quantile(matrix.smooth.selected$ratio, probs = 0.1)] <- 'k9'
matrix.smooth.selected$class[matrix.smooth.selected$ratio > quantile(matrix.smooth.selected$ratio, probs = 0.9)] <- 'k37'

table(matrix.smooth.selected$class)

smooth.E1.matrix <-left_join(smooth.E1.matrix,matrix.smooth.selected[, c('binID', 'class')],'binID' )

smooth.E1.matrix.diff <- smooth.E1.matrix[smooth.E1.matrix$class != 'ns',]


###### FINAL SCATTERPLOT FOR ENRICHED BINS
volcol <- c( '#f77979', '#7171f5',"#9C9B9B")
names(volcol) <- c("k37","k9","ns")

pdf(paste0('plots/scatter_smooth_enrich_', mark1, '_', mark2, '.pdf'))
ggplot(matrix.smooth.selected, aes(x=H3K37me3 , y=H3K9me3, color=class)) +
  geom_point(size = 1) +
  scale_colour_manual(values = volcol) + 
  ylim(-2,3) + xlim(-2,3)+ 
  xlab(mark1) + ylab(mark2) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 19))+
  geom_abline(slope = 1, intercept = 1, color = "black", linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1, color = "black", linetype = "dashed") + 
  geom_abline(slope = lm(matrix.smooth.selected[, mark1] ~ matrix.smooth.selected[, mark2])$coefficients[2],
              intercept = lm(matrix.smooth.selected[, mark1] ~ matrix.smooth.selected[, mark2])$coefficients[1]
              , color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.7) + 
  annotate("text", x = 0.1, y = 0.3, label = paste0('cor = ', round(cor(matrix.smooth.selected[, mark1],matrix.smooth.selected[, mark2]),3)), size = 5, color = "red") +
  ggtitle(paste0(mark1, ' vs. ', mark2, ' enriched in E1')) + 
  theme(plot.title = element_text(size = 18))

dev.off()



## BOXPLOTS FOR OMIC DATA IN ENRICHED BINS

pdf(file = paste0('plots/boxplot_smooth_enriched_', mark1, '_', mark2, '.pdf'), height = 10, width = 10)
colors <- c( '#f77979', '#7171f5',"#9C9B9B")
par(mfrow=c(3,4))
omics.names <- c("ATAC", "PolII", "Nascent RNA", "Oris",'H3K37me3','H3K9me3','H3K36me3', 'H3K37me1','H3K4me3','H3K27me3','H3K36me2','H3K9me1')
mins <- c(0,0,0,0, -0.1,-0.1,-0.1,-0.1, -0.1,-0.1,-0.1,-0.1)
maxs <- c(0.02,0.02,0.02,0.05, 0.3,0.3,0.05,0.05, 0.05, 0.05, 0.05, 0.05)

mins <- mins/2
maxs <- maxs/2

mins <- c(0,0,0,0, -2,-2,-1,-2, -1.5,-1.5,-3,-3)
maxs <- c(0.15,0.05,0.02,0.1,3, 3,1,1,1, 1, 1, 1, 1)

for (i in 1:length(omics.names)) {
  omic <- omics.names[i]
  boxplot(smooth.E1.matrix.diff[[omic]] ~ smooth.E1.matrix.diff$class, boxlwd = 2 , lwd = 2 ,  ylim=c(mins[i],maxs[i]), axes= F , xlab= "", ylab= "log2(CPM)", 
          col=colors, main = omic, las = 2, cex.axis = 0.7, outline = FALSE) 
  axis(2, cex.axis = 0.8)
  axis(1, at = 1:2, labels =c('K37 enrich','K9 enrich'), tick =F, las = 1, cex.axis = 0.9)
  text(x=1.5, y =  maxs[i] * 1, labels = paste0("p = ", wilcox.test(smooth.E1.matrix.diff[[omic]] ~ smooth.E1.matrix.diff$class)$p.value), cex = 0.8)
  #abline(h=0, lwd=2, lty=2)
}

dev.off()


