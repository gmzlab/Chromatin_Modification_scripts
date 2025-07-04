############################## LOADING AND PRE-PROCESSING COUNTS DATA ##############################

library(preprocessCore)
library(ggplot2)

# Loading raw counts table.
tab <- read.csv2('../results/chromHMM/counts/5kb_RawCount.txt' , header = T, sep = '\t', quote ="'", check.names = T, as.is = T)
colnames(tab) <- gsub('.bam', '', colnames(tab))
colnames(tab)[1] <- 'chr'
rownames(tab) <- paste('bin', 1:nrow(tab), sep='')

tab <- tab[grepl("^.{1,2}$", tab$chr), ]


# Getting only count data as numeric.
marks.bin.cnt <- tab[,-(1:3)]
for (i in 1:ncol(marks.bin.cnt)) { marks.bin.cnt[,i] <- as.numeric(marks.bin.cnt[,i]) }
marks.bin.cnt <- as.matrix(marks.bin.cnt)



# Calculating CPM data.
marks.bin.cpm <- t(t(marks.bin.cnt)/(colSums(marks.bin.cnt)/1e6))

# Getting positions data.
bin.bed <- tab[,1:3]
bin.bed$binID <- rownames(bin.bed)


# log2 transformation 

marks.bin.log2 <- log2(marks.bin.cpm +1)


# Input adjustment

marks.enrich <- marks.bin.log2[,1:14]

marks.enrich <- marks.enrich[,c('h3k37me3_B12_1','h3k9me3','h3k27me3','h3k9me2','h3k9me1','h3k36me1','h3k36me2','h3k36me3','h3k37me1_B12','h3k4me1','h3k4me2','h3k4me3','h3k9ac','h3k27ac')]

marks.enrich[,c("h3k27ac", "h3k27me3","h3k36me1","h3k36me2","h3k36me3","h3k4me1",
                "h3k4me2","h3k4me3","h3k9ac","h3k9me1","h3k9me2","h3k9me3")] <-  marks.bin.log2[,c("h3k27ac", "h3k27me3","h3k36me1","h3k36me2","h3k36me3","h3k4me1",
                                                                                                   "h3k4me2","h3k4me3","h3k9ac","h3k9me1","h3k9me2","h3k9me3")] - marks.bin.log2[, "input"]

marks.enrich[,"h3k37me1_B12"] <-  marks.bin.log2[,"h3k37me1_B12"] - marks.bin.log2[,"input1_B12"]
marks.enrich[,"h3k37me3_B12_1"] <-  marks.bin.log2[,"h3k37me3_B12_1"] - marks.bin.log2[,"input2_B12_1"]


colnames(marks.enrich) <-c('H3K37me3','H3K9me3','H3K27me3','H3K9me2','H3K9me1','H3K36me1','H3K36me2','H3K36me3','H3K37me1','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac')

boxplot(marks.enrich, ylim=c(-3,3), main = "Input Corrected", las = 2, cex.axis = 0.7, outline = FALSE)




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

# log2
omics.bin.log2 <- log2(omics.bin.cpm +1)


# INPUT CORRECTION


omics.enrich <- omics.bin.log2[,c("atac", "pol2", "h2az", "rna", "rna_nascent", "oris", "mcm3")]

omics.enrich[,"h2az"] <-  omics.bin.log2[,"h2az"] - omics.bin.log2[,"input2"]
omics.enrich[,"mcm3"] <-  omics.bin.log2[,"mcm3"] - omics.bin.log2[,"mcm3_input"]

colnames(omics.enrich) <- c("ATAC", "Pol II", "H2A.Z", "RNA", "Nascent RNA", "Oris", "MCM3")

boxplot(omics.enrich, ylim=c(-1,6), main = "Input corrected", las = 2, cex.axis = 0.7, outline = F)


summary(omics.enrich)



matrix <- left_join(marks.enrich, omics.enrich, 'binID')


############################ MARKS vs OMICS SCATTER #####


# Function for comparing selected marks with omics by scatterplots
plot_mark_vs_omic <- function(mark, omic) {
  # Initiate graph
  pdf(file = paste0('plots/scatterplots_', mark, '_vs_', omic, '.pdf'), height = 15, width = 10)
  par(mfrow = c(5, 1))
  
  create_plot <- function(data_mark, data_omic, title) {
    plot(data_mark[, mark], data_omic[, omic], xlim = c(-3, 3), ylim = c(-0.1, 6),
         main = paste(mark, 'vs.', omic, '~', title), xlab = mark, ylab = omic)
    
    correlation <- cor(data_mark[, mark], data_omic[, omic])
    text(x = 2.6, y = 5.5, labels = paste("Cor:", round(correlation, 2)), pos = 4)
    abline(lm(data_omic[, omic] ~ data_mark[, mark]), col = "red")
  }
  }
  
  # GENOME WIDE
  create_plot(marks.enrich, omics.enrich, "GENOME")
  
  # Activate states segmentation when aplicable
  ## E1
  #create_plot(marks.E1, omics.E1, "Heterochromatin")
  #
  ## E3
  #create_plot(marks.E3, omics.E3, "Facultative")
  #
  ## E4
  #create_plot(marks.E4, omics.E4, "Intermediate")
  #
  ## E5
  #create_plot(marks.E5, omics.E5, "Euchromatin")
  
  # Cerrar el dispositivo PDF
  dev.off()


mark <- "H3K37me3"
omic <- "atac"

plot_mark_vs_omic(mark,omic)

plot_mark_vs_omic("H3K37me3","pol2")
plot_mark_vs_omic("H3K37me3","rna_nascent")
plot_mark_vs_omic("H3K37me3","oris")


plot_mark_vs_omic("H3K9me3","pol2")
plot_mark_vs_omic("H3K9me3","rna_nascent")
plot_mark_vs_omic("H3K9me3","oris")





############################# CORRELATION BY QUANTILE #########

## INDIVIDUAL OMIC

plot_omic_quantiles <- function(mark, omic) {
  
# Adding quantile segregation as a factor to the main matrix
matrix[,paste(mark,".quantiles")] <- cut(matrix[,mark], breaks = quantile(matrix[,mark]), include.lowest = TRUE)
omic.quantiles <- split(matrix[,mark], matrix[,paste(mark,".quantiles")])

# Boxplot
boxplot(omic.quantiles, boxlwd = 2 , lwd = 2 , axes= F , ylim= c(-2, 5),xlab= "", ylab= "log2(CPM/Input)", 
        col=c('#97e4f0', '#19c6e0', '#3d8af5', '#021ef2' ), las = 2, cex.axis = 0.7,
        main = paste(omic, 'by', mark, 'quantiles'),  outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"), tick =F, las = 1, cex.axis = 1.3)

}

mark <- "H3K37me3"
omic <- "ATAC"

plot_omic_quantiles("H3K37me3", "a")


## ALL OMICS

plot_all_omics_quantiles <-  function(mark) {
  
  pdf(file = paste0('plots/boxplot_', mark, '_omics_quantile_test.pdf'), height = 15, width = 10)
  par(mfrow = c(3, 2))
  plot_omic_quantiles(mark, 'ATAC')
  plot_omic_quantiles(mark, 'Pol II')
  plot_omic_quantiles(mark, 'H2A.Z')
  plot_omic_quantiles(mark, 'RNA')
  plot_omic_quantiles(mark, 'Nascent RNA')
  plot_omic_quantiles(mark, 'MCM3')
  dev.off()
}

mark <- 'H3K37me3'

plot_all_omics_quantiles(mark)

plot_all_omics_quantiles('H3K37me3')
plot_all_omics_quantiles('H3K37me1')
plot_all_omics_quantiles('H3K9me3')
plot_all_omics_quantiles('H3K27me3')
plot_all_omics_quantiles('H3K36me3')
plot_all_omics_quantiles('H3K36me2')






############################# CORRELATION BY QUANTILE - NO 0 VALUES #####

# Removing null rows from the matrix to avoid interference with quantile segmentation

#matrix_no0 <- matrix[matrix[, mark] != 0, ]

matrix_no0 <- matrix[rowSums(matrix[,c(1:14,20:26)]) != 0, ]
  
nrow(matrix)-nrow(matrix_no0)

mark <- 'rna_nascent'

nrow(matrix[matrix[, mark] == 0, ])
nrow(matrix_no0[matrix_no0[, mark] == 0, ])


## Testing the results of null rows deletion

matrix_no0[,paste(mark,".quantiles")] <- cut(matrix_no0[,mark], breaks = quantile(matrix_no0[,mark]), include.lowest = TRUE)
omic.quantiles <- split(matrix_no0[,omic], matrix_no0[,paste(mark,".quantiles")])

# Boxplot
boxplot(omic.quantiles, boxlwd = 2 , lwd = 2 , axes= F , ylim= c(-2, 5),xlab= "", ylab= "log2(CPM/Input)", 
        col=c('#97e4f0', '#19c6e0', '#3d8af5', '#021ef2' ), las = 2, cex.axis = 0.7,
        main = paste(omic, 'by', mark, 'quantiles'),  outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"), tick =F, las = 1, cex.axis = 1.3)




## INDIVIDUAL OMICS

plot_omic_quantiles <- function(mark, omic, title, ylim = c(0, 5)) {
  
  # Adding quantile segregation as a factor to the main matrix
  matrix_no0[,paste(mark,".quantiles")] <- cut(matrix_no0[,mark], breaks = quantile(matrix_no0[,mark]), include.lowest = TRUE)
  omic.quantiles <- split(matrix_no0[,omic], matrix_no0[,paste(mark,".quantiles")])
  
  # Boxplot
  boxplot(omic.quantiles, boxlwd = 2 , lwd = 2 , axes= F , ylim= ylim,xlab= "", ylab= "log2(CPM+1)", 
          col=c('#97e4f0', '#19c6e0', '#3d8af5', '#021ef2' ), las = 2, cex.axis = 0.7,
          main = title, cex.main = 2 ,outline = FALSE) 
  axis(2, cex.axis = 0.8)
  axis(1, at = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"), tick =F, las = 1, cex.axis = 1.3)
  
}

mark <- "H3K37me3"
omic <- "oris"


cor <- cor(matrix[,mark],matrix[,omic])
print()
cor

plot_omic_quantiles(mark, omic)


## ALL OMICS

plot_all_omics_quantiles <-  function(mark) {
  
  pdf(file = paste0('plots/boxplot_', mark, '_omics_quantile.pdf'), height = 3.5, width = 14)
  par(mfrow = c(1, 4))
  plot_omic_quantiles(mark, 'ATAC', paste('ATAC'), ylim = c(0, 4))
  plot_omic_quantiles(mark, 'Pol.II', paste('Pol II'))
 # plot_omic_quantiles(mark, 'H2A.Z')
 # plot_omic_quantiles(mark, 'RNA')
  plot_omic_quantiles(mark, "Nascent.RNA", paste('Nascent RNA'))
  plot_omic_quantiles(mark, 'Oris', paste('Early Oris'))
  dev.off()
}

plot_all_omics_quantiles('H3K37me3')



# Select the desired marks for segregation
mark <- 'H3K37me3'
plot_all_omics_quantiles(mark)

plot_all_omics_quantiles('H3K37me1')
plot_all_omics_quantiles('H3K9me3')
plot_all_omics_quantiles('H3K27me3')
plot_all_omics_quantiles('H3K36me3')
plot_all_omics_quantiles('H3K36me2')




matrix_no0[,paste(mark,".quantiles")] <- cut(matrix_no0[,mark], breaks = quantile(matrix_no0[,mark]), include.lowest = TRUE)
omic.quantiles <- split(matrix_no0[,"Oris"], matrix_no0[,paste(mark,".quantiles")])

# Boxplot
boxplot(omic.quantiles, boxlwd = 2 , lwd = 2 , axes= F , ylim= c(-2, 5),xlab= "", ylab= "log2(CPM/Input)", 
        col=c('#97e4f0', '#19c6e0', '#3d8af5', '#021ef2' ), las = 2, cex.axis = 0.7,
        main = paste(omic, 'by', mark, 'quantiles'),  outline = FALSE) 
axis(2, cex.axis = 0.8)
axis(1, at = 1:4, labels = c("Q1", "Q2", "Q3", "Q4"), tick =F, las = 1, cex.axis = 1.3)



