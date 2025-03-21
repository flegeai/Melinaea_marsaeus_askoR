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
parameters$contrast_file = "Contrasts.txt"            # matrix of different contrasts desired
parameters$sample_file = "Samples.csv"   # file describing the samples
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
resDEG<-DEanalysis(asko_norm, data, asko_data, parameters)

###### UpsetR Graphs #####
# My list
parameters$upset_list = c("Mmp_1vsMmr_1-Mmp_2vsMmr_2-Mmp_aavsMmr_aa",
                          "WD1_MmrvsWD2_Mmr-WD1_MmpvsWD2_Mmp-WD1vsWD2")

# graphs type "all"
parameters$upset_basic = "all"
parameters$upset_type = "all"
UpSetGraph(resDEG, data, parameters)

# graphs type "mixed"
parameters$upset_basic = "mixed"
parameters$upset_type = "mixed"
UpSetGraph(resDEG, data, parameters)

# graphs type "up"
parameters$upset_basic = "up"
parameters$upset_type = "up"
UpSetGraph(resDEG, data, parameters)

# graphs type "down"
parameters$upset_basic = "down"
parameters$upset_type = "down"
UpSetGraph(resDEG, data, parameters)

##### Enrichment Analysis #####
# Parameters
#----------------------------------------------------------------------


parameters$GO_threshold = 0.05               # the significant threshold used to filter p-values (default 0.05)
parameters$GO_max_top_terms = 10             # the maximum number of GO terms plot (default 10)
parameters$GO_min_num_genes = 10             # the minimum number of genes for each GO terms (default 10)
parameters$GO_min_sig_genes = 0              # the minimum number of significant gene(s) behind the enriched GO-term (default 0)
parameters$GO = "both"                       # gene set chosen for analysis 'up', 'down', 'both', or NULL (default NULL)
parameters$GO_algo = "weight01"              # algorithms for runTest function ("classic", "elim", "weight", "weight01", "lea", "parentchild") (default weight01)
parameters$GO_stats = "fisher"               # statistical tests for runTest function ("fisher", "ks", "t", "globaltest", "sum", "ks.ties") (default fisher)
parameters$Ratio_threshold = 1               # the min ratio for display GO in graph (default 0)
GOenrichment(resDEG, data, parameters)
