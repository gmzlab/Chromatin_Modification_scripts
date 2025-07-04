

######################################## LOADING RAW COUNTS DATA #############


# Original merged .bams from 14 marks + 3 inputs, quantified and compiled with deeptools multiBamSummary
tab <- read.csv2('../results/chromHMM/counts/5kb_rawCount.txt' , header = T, sep = '\t', quote ="'", check.names = T, as.is = T)
colnames(tab) <- gsub('.bam', '', colnames(tab))
colnames(tab)[1] <- 'chr'
rownames(tab) <- paste('bin', 1:nrow(tab), sep='')

# Removing extra scaffolds
tab <- tab[grepl("^.{1,2}$", tab$chr), ]
tab <- tab[tab$chr != 'MT',]

table(tab$chr)

# 2 + 2 replicates of 37me1 B10/B12 + 4 respective inputs for data transformation comparison
tab2 <- read.csv2('../results/chromHMM/counts/replicates_5kb_rawCount.txt' , header = T, sep = '\t', quote ="'", check.names = T, as.is = T)
colnames(tab2)[1] <- 'chr'
rownames(tab2) <- paste('bin', 1:nrow(tab2), sep='')

tab2 <- tab2[grepl("^.{1,2}$", tab2$chr), ]
tab2 <- tab2[tab2$chr != 'MT',]

table(tab2$chr)


# Getting positions data.
bin.bed <- tab[,1:3]
bin.bed$binID <- rownames(bin.bed)


# Combining count data from both tables as numeric

marks.bin.cnt <- cbind(tab[,-(1:3)],tab2[,-(1:3)])

for (i in 1:ncol(marks.bin.cnt)) { marks.bin.cnt[,i] <- as.numeric(marks.bin.cnt[,i]) }
marks.bin.cnt <- as.matrix(marks.bin.cnt)


## Selecting and reordering columns

# Only merged 
# marks.bin.cnt <- marks.bin.cnt[,c('h3k37me3_B12_1','h3k9me3','h3k27me3','h3k9me2','h3k9me1','h3k36me1','h3k36me2','h3k36me3','h3k37me1_B12','h3k4me1','h3k4me2','h3k4me3','h3k9ac','h3k27ac')]

# Including Replicates comparison
marks.bin.cnt <- marks.bin.cnt[,c('h3k37me3_B12_1','h3k9me3','h3k27me3','h3k9me2','h3k9me1','h3k36me1','h3k36me2',
                                  'h3k36me3','h3k4me1','h3k4me2','h3k4me3','h3k9ac','h3k27ac', 
                                  'h3k37me1_B12', 'h3k37me1_B10_1', 'h3k37me1_B12_1', 'h3k37me1_B10_2', 'h3k37me1_B12_2',
                                  'input1_B12', 'input2_B12_1', 'input', 'input1_B10_1', 'input1_B12_1', 'input1_B10_2', 'input1_B12_2')]


colnames(marks.bin.cnt) <-c('H3K37me3','H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2',
                            'H3K36me3','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac',
                            'H3K37me1 merge','37me1 B10 1','37me1 B12 1','37me1 B10 2','37me1 B12 2',
                            'input 37me1 merge', 'input 37me3', 'input rest', 'input 37me1 B10 1','input 37me1 B12 1','input 37me1 B10 2','input 37me1 B12 2')


############################ PRE-PROCESSING ################


# 1. CALCULATING CPM
marks.bin.cpm <- t(t(marks.bin.cnt)/(colSums(marks.bin.cnt)/1e6))

boxplot(marks.bin.cpm, ylim=c(0,7), main = "cpm", las = 2, cex.axis = 0.7, outline = FALSE)


# 2. REMOVING NULL ROWS
marks.bin.cpm.no0 <- marks.bin.cpm[rowSums(marks.bin.cpm) != 0, ]

boxplot(marks.bin.cpm.no0, ylim=c(0,7), main = "no0", las = 2, cex.axis = 0.7, outline = FALSE)


# 3. Log2 + 1
marks.bin.log2 <- log2(marks.bin.cpm.no0 +1)

boxplot(marks.bin.log2, ylim=c(0,7), main = "log2", las = 2, cex.axis = 0.7, outline = FALSE)


# 4. INPUT CORRECTION
# Only merged samples
marks.corrected <- marks.bin.log2[,1:14]

#Including replicates for comparison
#marks.corrected <- marks.bin.log2[,1:25]

# GSE175752 dataset
marks.corrected[,c('H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2',
                  'H3K36me3','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')] <-  marks.bin.log2 [,c('H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2',
                                                                                                       'H3K36me3','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')] - marks.bin.log2[,'input rest']
# B10/B12 ChIP-seq 
marks.corrected[,'H3K37me1 merge'] <- marks.bin.log2[,'H3K37me1 merge'] - marks.bin.log2[,'input 37me1 merge']
marks.corrected[,'H3K37me3'] <- marks.bin.log2[,'H3K37me3'] - marks.bin.log2[,'input 37me3']

marks.corrected[,'37me1 B10 1'] <- marks.bin.log2[,'37me1 B10 1'] - marks.bin.log2[,'input 37me1 B10 1']
marks.corrected[,'37me1 B12 1'] <- marks.bin.log2[,'37me1 B12 1'] - marks.bin.log2[,'input 37me1 B12 1']
marks.corrected[,'37me1 B10 2'] <- marks.bin.log2[,'37me1 B10 2'] - marks.bin.log2[,'input 37me1 B10 2']
marks.corrected[,'37me1 B12 2'] <- marks.bin.log2[,'37me1 B12 2'] - marks.bin.log2[,'input 37me1 B12 2']


boxplot(marks.corrected, ylim=c(-3,3), main = "input corrected", las = 2, cex.axis = 0.7, outline = FALSE)


# Reordering again

marks.corrected <- marks.corrected[,c('H3K37me3','H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2',
                                      'H3K36me3', 'H3K37me1 merge', 'H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')]


colnames(marks.corrected) <-c('H3K37me3','H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2',
                            'H3K36me3', 'H3K37me1', 'H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')




# 4.5 INPUT CORRECTED (no log2)
# Testing for alternative without log2 transformation, so input counts are divided instead of substracted
# Pseoudocount transformation is added to avoid dividing by 0
#pseudocount <- 1
#
#marks.corrected.nolog <- marks.bin.cpm.no0[,1:18] 
#marks.bin.cpm.no0.pseudo <- marks.bin.cpm.no0 + pseudocount
#
#
#marks.corrected.nolog[,c('H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2',
#                   'H3K36me3','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')] <-  marks.bin.cpm.no0.pseudo [,c('H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2',
#                                                                                                        'H3K36me3','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')] / marks.bin.cpm.no0.pseudo[,'input rest']
#marks.corrected.nolog[,'H3K37me1 merge'] <- marks.bin.cpm.no0.pseudo[,'H3K37me1 merge'] / marks.bin.cpm.no0.pseudo[,'input 37me1 merge']
#marks.corrected.nolog[,'H3K37me3'] <- marks.bin.cpm.no0.pseudo[,'H3K37me3'] / marks.bin.cpm.no0.pseudo[,'input 37me3']
#
#marks.corrected.nolog[,'37me1 B10 1'] <- marks.bin.cpm.no0.pseudo[,'37me1 B10 1'] / marks.bin.cpm.no0.pseudo[,'input 37me1 B10 1']
#marks.corrected.nolog[,'37me1 B12 1'] <- marks.bin.cpm.no0.pseudo[,'37me1 B12 1'] / marks.bin.cpm.no0.pseudo[,'input 37me1 B12 1']
#marks.corrected.nolog[,'37me1 B10 2'] <- marks.bin.cpm.no0.pseudo[,'37me1 B10 2'] / marks.bin.cpm.no0.pseudo[,'input 37me1 B10 2']
#marks.corrected.nolog[,'37me1 B12 2'] <- marks.bin.cpm.no0.pseudo[,'37me1 B12 2'] / marks.bin.cpm.no0.pseudo[,'input 37me1 B12 2']
#
#
#boxplot(marks.corrected.nolog, ylim=c(-1,4), main = "input corrected no log2", las = 2, cex.axis = 0.7, outline = FALSE)



########################### PCA #############################

# Exploratory PCA not suitable for final figures
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(ggrepel))

pca <- PCA(t(marks.bin.cpm), graph = FALSE)
fviz_pca_ind(pca, pointsize=2, pointshape=21,fill="black",
             repel = TRUE, show_legend=TRUE,show_guide=FALSE, title = 'CPM')

pca <- PCA(t(marks.bin.cpm.no0), graph = FALSE)
fviz_pca_ind(pca, pointsize=2, pointshape=21,fill="black",
             repel = TRUE, show_legend=TRUE,show_guide=FALSE, title = 'no0')

pca <- PCA(t(marks.bin.log2), graph = FALSE)
fviz_pca_ind(pca, pointsize=2, pointshape=21,fill="black",
             repel = TRUE, show_legend=TRUE,show_guide=FALSE, title = 'log2')

pca <- PCA(t(marks.corrected), graph = FALSE)
fviz_pca_ind(pca, pointsize=2, pointshape=21,fill="black",
             repel = TRUE, show_legend=TRUE,show_guide=FALSE, title = 'input corrected')

pca <- PCA(t(marks.corrected.nolog), graph = FALSE)
fviz_pca_ind(pca, pointsize=2, pointshape=21,fill="black",
             repel = TRUE, show_legend=TRUE,show_guide=FALSE, title = 'input corrected no log 2')




############################ HEATMAP CORRELATION CLUSTERING ################
suppressMessages(library(heatmaply))
suppressMessages(library(ggplot2))

# Use the desired matrix for clustering
cor.matrix <- cor(marks.bin.cpm)
cor.matrix <- cor(marks.bin.cpm.no0)
cor.matrix <- cor(marks.bin.log2)
cor.matrix <- cor(marks.corrected)
cor.matrix <- cor(marks.corrected.nolog)


# Viewer + Dendogram
heatmaply_cor(cor.matrix, show_dendrogram = TRUE, main = "Histone PTM Correlation")



# With explicit values
ggplot(reshape2::melt(cor.matrix), aes(Var2, Var1, fill = value))+ggtitle("input corrected")+ geom_tile(color = 'white')+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") + 
  geom_text(aes(Var2, Var1, label = round(value, 2)), size = 2.2) +theme_minimal() +
  coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1)) +
  labs(x = "", y = "")

