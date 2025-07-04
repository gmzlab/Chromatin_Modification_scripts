

######################################## LOADING RAW COUNTS DATA ######


# Loading counting matrix whith 2 replicates for each mark

tab <- read.csv2('../Data/WT/5kb_replicates_RawCount.txt' , header = T, sep = '\t', quote ="'", check.names = T, as.is = T)

colnames(tab) <- gsub('.bam', '', colnames(tab))
colnames(tab)[1] <- 'chr'
rownames(tab) <- paste('bin', 1:nrow(tab), sep='')

# Removing extra chromosomes
tab <- tab[grepl("^.{1,2}$", tab$chr), ]
tab <- tab[tab$chr != 'MT',]

table(tab$chr)

# Getting positions data.
bin.bed <- tab[,1:3]
bin.bed$binID <- rownames(bin.bed)


# Combining count data as numeric

marks.bin.cnt <- cbind(tab[,-(1:3)])

for (i in 1:ncol(marks.bin.cnt)) { marks.bin.cnt[,i] <- as.numeric(marks.bin.cnt[,i]) }
marks.bin.cnt <- as.matrix(marks.bin.cnt)


## Reordering columns for replicates comparison
marks.bin.cnt <- marks.bin.cnt[,c('h3k37me3_B12_1',
                                  'h3k9me3_1', 'h3k9me3_2',
                                  'h3k27me3_1', 'h3k27me3_2', 
                                  'h3k9me2_1', 'h3k9me2_2',
                                  'h3k9me1_1', 'h3k9me1_2',
                                  'h3k36me1_1', 'h3k36me1_2',
                                  'h3k36me2_1', 'h3k36me2_2',
                                  'h3k36me3_1', 'h3k36me3_2',
                                  'h3k4me1_1', 'h3k4me1_2',
                                  'h3k4me2_1', 'h3k4me2_2',
                                  'h3k4me3_1', 'h3k4me3_2',
                                  'h3k9ac_1', 'h3k9ac_2',
                                  'h3k27ac_1', 'h3k27ac_2',
                                  'h3k37me1_B12_1', 'h3k37me1_B12_2',
                                  'input_1', 'input_2',
                                  'input1_B12_1', 'input1_B12_2',
                                  'input2_B12_1')]


colnames(marks.bin.cnt) <-c('H3K37me3 1',
                            'H3K9me3 1', 'H3K9me3 2',
                            'H3K27me3 1', 'H3K27me3 2',
                            'H3K9me2 1', 'H3K9me2 2',
                            'H3K9me1 1', 'H3K9me1 2',
                            'H3K36me1 1', 'H3K36me1 2',
                            'H3K36me2 1', 'H3K36me2 2',
                            'H3K36me3 1', 'H3K36me3 2',
                            'H3K4me1 1', 'H3K4me1 2',
                            'H3K4me2 1', 'H3K4me2 2',
                            'H3K4me3 1', 'H3K4me3 2',
                            'H3K9ac 1', 'H3K9ac 2',
                            'H3K27ac 1', 'H3K27ac 2',
                            'H3K37me1 1', 'H3K37me1 2',
                            'input_1', 'input_2',
                            'input1_B12_1', 'input1_B12_2',
                            'input2_B12_1')

############################ DATA PRE-PROCESSING ################


# 1. CALCULATING CPM
marks.bin.cpm <- t(t(marks.bin.cnt)/(colSums(marks.bin.cnt)/1e6))

# 2. REMOVING NULL ROWS
marks.bin.cpm.no0 <- marks.bin.cpm[rowSums(marks.bin.cpm) != 0, ]

# 3. Log2 + 1
marks.bin.log2 <- log2(marks.bin.cpm.no0 +1)

# 4. INPUT CORRECTION
marks.corrected <- marks.bin.log2[,1:27]

marks.corrected[,c('H3K9me3 1',
                   'H3K27me3 1',
                   'H3K9me2 1',
                   'H3K9me1 1',
                   'H3K36me1 1',
                   'H3K36me2 1',
                   'H3K36me3 1',
                   'H3K4me1 1',
                   'H3K4me2 1',
                   'H3K4me3 1',
                   'H3K9ac 1',
                   'H3K27ac 1')] <-  marks.bin.log2 [,c('H3K9me3 1',
                                                         'H3K27me3 1',
                                                         'H3K9me2 1',
                                                         'H3K9me1 1',
                                                         'H3K36me1 1',
                                                         'H3K36me2 1',
                                                         'H3K36me3 1',
                                                         'H3K4me1 1',
                                                         'H3K4me2 1',
                                                         'H3K4me3 1',
                                                         'H3K9ac 1',
                                                         'H3K27ac 1')] - marks.bin.log2[,'input_1']

marks.corrected[,c('H3K9me3 2',
                   'H3K27me3 2',
                   'H3K9me2 2',
                   'H3K9me1 2',
                   'H3K36me1 2',
                   'H3K36me2 2',
                   'H3K36me3 2',
                   'H3K4me1 2',
                   'H3K4me2 2',
                   'H3K4me3 2',
                   'H3K9ac 2',
                   'H3K27ac 2')] <-  marks.bin.log2 [,c('H3K9me3 2',
                                                        'H3K27me3 2',
                                                        'H3K9me2 2',
                                                        'H3K9me1 2',
                                                        'H3K36me1 2',
                                                        'H3K36me2 2',
                                                        'H3K36me3 2',
                                                        'H3K4me1 2',
                                                        'H3K4me2 2',
                                                        'H3K4me3 2',
                                                        'H3K9ac 2',
                                                        'H3K27ac 2')] - marks.bin.log2[,'input_2']



marks.corrected[,'H3K37me1 1'] <- marks.bin.log2[,'H3K37me1 1'] - marks.bin.log2[,'input1_B12_1']
marks.corrected[,'H3K37me1 2'] <- marks.bin.log2[,'H3K37me1 2'] - marks.bin.log2[,'input1_B12_2']

marks.corrected[,'H3K37me3 1'] <- marks.bin.log2[,'H3K37me3 1'] - marks.bin.log2[,'input2_B12_1']


boxplot(marks.corrected, ylim=c(-3,3), main = "input corrected replicates", las = 2, cex.axis = 0.7, outline = FALSE)


########################### PCA #############################

suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(ggrepel))
suppressMessages(library(dplyr))


pca <- PCA(t(marks.corrected), graph = FALSE)


# Defining the marks as ordered factors to correctly tag the points in the plot
experimental.design <- data.frame(sample = colnames(marks.corrected), condition = colnames(marks.corrected))

experimental.design$condition <- sub(" \\d+$", "", experimental.design$condition)


label.positions <- !duplicated(experimental.design$condition)
labels <- ifelse(label.positions, experimental.design$condition, "")

experimental.design$condition <- factor(experimental.design$condition,
                                        levels = unique(experimental.design$condition))                                       

# PCA ~ Tagging all replicates
fviz_pca_ind(pca, col.ind = experimental.design$condition, 
             pointsize=2, pointshape=21, fill.ind = experimental.design$condition,
             repel = TRUE,
             legend.title="Histone Marks",
             title= "Principal Components Analysis",
             show_legend= TRUE, show_guide=TRUE,
             xlab = paste0('PC1 (',round(pca$eig[1,2],1),'%)'),
             ylab = paste0('PC2 (',round(pca$eig[2,2],1),'%)')) 



# Tagging only conditions
pca_coords <- as.data.frame(pca$ind$coord)
pca_coords$condition <- experimental.design$condition
pca_coords$condition <- factor(pca_coords$condition, levels = unique(experimental.design$condition))

# Obtaining centroids
centroids <- pca_coords %>%
  group_by(condition) %>%
  summarise(Dim.1 = mean(Dim.1),
            Dim.2 = mean(Dim.2))

# Assign condition's order in centroids
centroids$condition <- factor(centroids$condition, levels = levels(pca_coords$condition))

# PCA ~  Tagging centroids only, no legend 
pdf(file = paste0('plots/pca_replicates.pdf'), height = 5, width = 7)

fviz_pca_ind(pca,
             col.ind = experimental.design$condition,
             fill.ind = experimental.design$condition,
             pointsize = 2, pointshape = 21,
             label = "none",
             legend.title = "Histone Marks",
             title = "Principal Components Analysis",
             show_legend = T, show_guide = F,
             xlab = paste0('PC1 (', round(pca$eig[1,2],1), '%)'),
             ylab = paste0('PC2 (', round(pca$eig[2,2],1), '%)')) + 
            geom_text_repel(data = centroids,
              aes(x = Dim.1, y = Dim.2, label = condition),
              size = 4,
              box.padding = 0.4,
              show.legend = FALSE) + 
            theme(legend.position = "none")

dev.off()
