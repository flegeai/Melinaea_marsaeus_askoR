setwd("~/tmp/Melinaea_marsaeus_askoR/ScriptR")
rm(list=ls())
# source AskoR.R file
source("AskoR.R")
parameters<-Asko_start()

##############################################
##                Parameters                ##
##############################################

# Data and input files descriptions
#--------------------------------------------------------------------------
# WARNING: All the input files must be in the same folder
#          called "input" (case sensitive)!
#--------------------------------------------------------------------------
parameters$dir_path = "../data"
parameters$analysis_name = "DE"             # output directory name (default DE_analysis, do not put space!)
#parameters$fileofcount = "CountsMatrix.txt"           # matrix of count for all samples/conditions
#parameters$sep = "\t"                                 # field separator for count files or count matrix
#parameters$annotation = "Genes_annotations.txt"       # file containing the functional annotations of each gene
parameters$geneID2GO_file = "MmEP.gos"      # GO annotation files
parameters$contrast_file = "Contrasts_noAA.txt"            # matrix of different contrasts desired
parameters$sample_file = "Samples_noAA.csv"   # file describing the samples
parameters$col_genes=1
parameters$col_counts=3
parameters$rm_sample = c("Mmr37_w1")             # bad sample(s) !

# Options for data processing and their analyzes
#--------------------------------------------------------------------------
parameters$threshold_cpm = 0.5                       # CPM's threshold (default 0.5)
parameters$replicate_cpm = 3                          # Minimum number of replicates (default 3)
parameters$threshold_FDR = 0.05                       # FDR threshold (default 0.05)
parameters$threshold_logFC = 0                        # logFC threshold (default 1)
parameters$normal_method = "TMM"                      # normalization method (TMM/RLE/upperquartile/none) (default TMN)
parameters$p_adj_method = "BH"                        # p-value adjust method (holm/hochberg/hommel/bonferroni/BH/BY/fdr/none) (default fdr)
parameters$glm = "qlf"                                # GLM method (lrt/qlf) (default qlf)
parameters$logFC = TRUE                               # logFC in the summary table (default TRUE)
parameters$FC = TRUE                                  # FC in the summary table (default TRUE)
parameters$logCPM = FALSE                             # logCPm in the summary table (default FALSE)
parameters$FDR = TRUE                                 # FDR in the summary table (default TRUE)
parameters$LR = FALSE                                 # LR in the summary table (default FALSE)
parameters$Sign = TRUE                                # Significance (1/0/-1) in the summary table (default TRUE)
parameters$Expression = TRUE                          # Significance expression in the summary table (default TRUE)
parameters$mean_counts = FALSE                        # Mean counts in the summary table (default TRUE)
parameters$norm_counts = FALSE                        # Generate files with mormalized counts

# for legend of density plot
#-----------------------------------
parameters$densbotmar = 15                            # Set bottom margin of density plot to help position the legend (default 20)
parameters$densinset = 0.20                           # Set position the legend in bottom density graphe (default 0.45)
parameters$legendcol = 6                              # Set numbers of column for legends (default 6)

# Visualization of results from differential expression analyses
#-----------------------------------
parameters$plotMD = TRUE                              # Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
parameters$plotVO = TRUE                              # Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)
parameters$glimMD = TRUE                              # Glimma - Interactif Mean-Difference Plot of Expression Data (aka MA plot) (default FALSE)
parameters$glimVO = TRUE                              # Glimma - Interactif Volcano plot for a specified coefficient/contrast of a linear model (default FALSE)

########################################
##  Loading the data from the samples ##
########################################
##### load data #####
cat("\nRun DE analysis\n")
data<-loadData(parameters)
cat("\n\nChecking Data content:\n")
data$samples
data$contrast
data$design

head(data$dge$counts,n=4)

cat("Total number of genes : ", dim(data$dge$counts)[1], "\n")
cat("Total number of samples : ", dim(data$dge$counts)[2], "\n\n")
cat("summary of CPM by samples\n")
summary(edgeR::cpm(data$dge))
cat("\n")

##### asko files #####
asko_data<-asko3c(data, parameters)
cat("\nChecking Asko Data : condition, contrast, context.\n")
asko_data$condition ; cat("\n")
asko_data$contrast  ; cat("\n")
asko_data$context   ; cat("\n")


write.table(edgeR::cpm(data$dge), file="CPM_all.txt")

##### filtering #####
cat("\nFiltering genes with more than ", parameters$threshold_cpm, " CPM in ",parameters$replicate_cpm,"samples\n")
asko_filt<-GEfilt(data, parameters)
cat("Total number of filtered genes : ", dim(asko_filt$counts)[1], "\n\n")

##### normalization #####
asko_norm<-GEnorm(asko_filt, asko_data, data, parameters)

##### correlation #####
GEcorr(asko_norm,parameters)

##### DGE analysis #####
cat("\n\nDifferential expressions analysis\n")
species<-factor(data$dge$samples$species)
stage<-factor(data$dge$samples$stage)
design <- model.matrix(~species + stage + species:stage, data=asko_norm)
design
normGEdisp <- estimateGLMCommonDisp(asko_norm, design)
normGEdisp <- estimateGLMTrendedDisp(normGEdisp, design)
normGEdisp <- estimateGLMTagwiseDisp(normGEdisp, design)
fit <- glmQLFit(normGEdisp, design, robust = TRUE)
test2<-glmQLFTest(fit, coef=2)
rileyi=decideTestsDGE(test2, adjust.method = 'BH', lfc=0, p.value=0.05)
write.csv(rileyi, "rileyi.csv")
summary(rileyi)
test3<-glmQLFTest(fit, coef=3)
stage=decideTestsDGE(test3, adjust.method = 'BH', lfc=0, p.value=0.05)
write.csv(stage, "stage.csv")
summary(stage)
test4<-glmQLFTest(fit, coef=4)
interaction=decideTestsDGE(test4, adjust.method = 'BH', lfc=0, p.value=0.05)
write.csv(interaction, "interaction.csv")
summary(interaction)
