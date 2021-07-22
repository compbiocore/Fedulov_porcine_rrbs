#######################
# title: "RRBS"
# author: "Jordan"
# date: "7/15/2021"
# output: html_document
#######################


# Note: Do liberal, moderate, and conservative approaches - 3 analyses with results 


# Load necessary packages (make sure they are installed) 
library(biomaRt)
#library(DSS)
library(bsseq)
library(dplyr)
#library(AnnotationDbi)
library(GenomicRanges)
#library(clusterProfiler)
#library(methyAnalysis)
#library(DESeq2)
library(edgeR)
#library(PCAtools)
library(openxlsx)

# Notes: 
# We ran a negative binomial GLM model with a one-way layout, i.e., treated each combination of treatments as groups and performed contrasts with these groups; this fit the research questions asked 


# Bring in targets file which lays out design and analytic setup 
getwd()
setwd("/Users/jordan/Desktop/rrbs_r")
library(readxl)
targets_2 <- read_excel("target.xlsx")
names(targets_2)
targets_2<-as.data.frame(targets_2[,-c(3,4,8)])
head(targets_2)
# Reference target file 
# targets_jos <- read.table('targets_2.txt', sep = '\t', header = T) #targets and targets_2 change additive and control factors

# import those files using readBismark2DGE (can be .cov or .cov.gz)
files <- targets_2$file #these are the names of the cov files 
setwd("/Users/jordan/Desktop/rrbs_r/cov_files")
yall <- readBismark2DGE(files, sample.names=targets_2$sample)

# filter to include only some chromosomes and sort them (we skip this for now)

# Question for Joselynn, should we remove mitochondrial genes "MT"? This is usually done according to Chen et al. 

# NOTE: Remove everything but 1:18 and X and Y 

# Look at how many - if couple hundred filter out! 

# yall <- yall[yall$genes[,'Chr'] == 1 | yall$genes[,'Chr'] == 2| yall$genes[,'Chr'] == 3| yall$genes[,'Chr'] == 4| yall$genes[,'Chr'] == 5| yall$genes[,'Chr'] == 6| yall$genes[,'Chr'] == 7| yall$genes[,'Chr'] == 8| yall$genes[,'Chr'] == 9| yall$genes[,'Chr'] == 10| yall$genes[,'Chr'] == 11| yall$genes[,'Chr'] == 12| yall$genes[,'Chr'] == 13| yall$genes[,'Chr'] == 14| yall$genes[,'Chr'] == 15| yall$genes[,'Chr'] == 16| yall$genes[,'Chr'] == 17| yall$genes[,'Chr'] == 18| yall$genes[,'Chr'] == 19| yall$genes[,'Chr'] == 'X',]
# Sort and re-organize genes data 
backup<-yall
sort_order <- c(1:18, 'X', 'Y'); sort_order
library(Hmisc)
values <- yall$genes[which(yall$genes$Chr %nin% sort_order),] 
values <- as.character(values$Chr) 
values <- unique(values); values # unique chromosomes that are not 1 through 18 and X and Y 
library(stringr)
values <- str_sort(values); values
sort_order <- c(1:18, 'X', 'Y', values); sort_order
yall$genes$Chr <- factor(yall$genes$Chr, levels=sort_order)
o <- order(yall$genes$Chr, yall$genes$Locus)
yall <- yall[o,] # reordering rows 


#####################################################################################################

# Right now, we did NOT annotate the CpG loci with the identity of the nearest gene! 

#add info about nearest TSS - going to have to figure this out later and tweak this for S_scrofa 

#BiocManager::install("org.Ss.eg.db")
# library(org.Ss.eg.db)
# right now you have yall$genes vars: Chr = Chromosome and Locus 
# below adds ID, Symbol, Strand, Distance, Width 
# library(AnnotationHub)
# ah <- AnnotationHub()
# orgdb <- query(ah, c("Sus scrofa", "OrgDb")) 
# TSS <- nearestTSS(yall$genes$Chr, yall$genes$Locus, species="Mm") # need org.Ss.egCHRLOC object 
# yall$genes$EntrezID <- TSS$gene_id # gene entrez id
# yall$genes$Symbol <- TSS$symbol # gene symbol
# yall$genes$Strand <- TSS$strand
# yall$genes$Distance <- TSS$distance # integer vector giving distance to nearest TSS. Positive values means that the TSS is downstream of the locus, negative values means that it is upstream. Gene body loci will therefore have negative distances and promotor loci will have positive.
# yall$genes$Width <- TSS$width # genomic width of the gene

########################################################################################################

# Gene filtering and normalization 

# make factor levels of Me and Un repeated along the length of the columns 
Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
Me <- yall$counts[, Methylation=="Me"]
Un <- yall$counts[, Methylation=="Un"]

# combine methylated and unmethylated reads to get total coverage
Coverage <- Me + Un

# Get a list of CpGs with at least 60 counts (methylated and unmethylated) in at least 3 samples - too strict
#HasCoverage <- rowSums(Coverage >= 60) == 3 # requires that sum across methylated and unmethylated counts is at least 60 across three samples 

# Get a list of CpGs with at least 8 counts (methylated and unmethylated) in all samples - conservative rule of thumb proposed by Chen et al. - too conservative 
HasCoverage <- rowSums(Coverage >= 8) == 18 # requires that sum across methylated and unmethylated counts is at least 60 across three samples 

# List of CpGs that are not 0 both
# this filters out CpG that are never methylated or always methylated as this provides no info about differential methylation 
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0 

# Check sample size after filtering requirements are met 
table(HasCoverage, HasBoth)

# apply those coverage filters
y <- yall[HasCoverage & HasBoth,, keep.lib.sizes=FALSE] # keeps rows that don't have all zeros and have at least 3 Me and Un combinations that total at least 60 - this drops a lot of data! 

#assign library sizes to each sample to be a sum of methylated and unmethylated counts
TotalLibSize <- y$samples$lib.size[Methylation=="Me"] + y$samples$lib.size[Methylation=="Un"]
y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples

# make some output tables for later..
counts_and_genes <- inner_join(tibble::rownames_to_column(as.data.frame(y$counts), var = 'loci'), tibble::rownames_to_column(as.data.frame(y$genes), var = 'loci'))
counts_and_genes$methylation_proportion_1 <- (counts_and_genes$'RH-1-Me') / (counts_and_genes$'RH-1-Me' + counts_and_genes$'RH-1-Un')
counts_and_genes$methylation_proportion_2 <- (counts_and_genes$'RH-2-Me') / (counts_and_genes$'RH-2-Me' + counts_and_genes$'RH-2-Un')
counts_and_genes$methylation_proportion_3 <- (counts_and_genes$'RH-3-Me') / (counts_and_genes$'RH-3-Me' + counts_and_genes$'RH-3-Un')
counts_and_genes$methylation_proportion_4 <- (counts_and_genes$'RH-4-Me') / (counts_and_genes$'RH-4-Me' + counts_and_genes$'RH-4-Un')
counts_and_genes$methylation_proportion_5 <- (counts_and_genes$'RH-5-Me') / (counts_and_genes$'RH-5-Me' + counts_and_genes$'RH-5-Un')
counts_and_genes$methylation_proportion_6 <- (counts_and_genes$'RH-6-Me') / (counts_and_genes$'RH-6-Me' + counts_and_genes$'RH-6-Un')
counts_and_genes$methylation_proportion_7 <- (counts_and_genes$'RH-7-Me') / (counts_and_genes$'RH-7-Me' + counts_and_genes$'RH-7-Un')
counts_and_genes$methylation_proportion_8 <- (counts_and_genes$'RH-8-Me') / (counts_and_genes$'RH-8-Me' + counts_and_genes$'RH-8-Un')
counts_and_genes$methylation_proportion_10 <- (counts_and_genes$'RH-10-Me') / (counts_and_genes$'RH-10-Me' + counts_and_genes$'RH-10-Un')
counts_and_genes$methylation_proportion_11 <- (counts_and_genes$'RH-11-Me') / (counts_and_genes$'RH-11-Me' + counts_and_genes$'RH-11-Un')
counts_and_genes$methylation_proportion_12 <- (counts_and_genes$'RH-12-Me') / (counts_and_genes$'RH-12-Me' + counts_and_genes$'RH-12-Un')
counts_and_genes$methylation_proportion_14 <- (counts_and_genes$'RH-14-Me') / (counts_and_genes$'RH-14-Me' + counts_and_genes$'RH-14-Un')
counts_and_genes$methylation_proportion_15 <- (counts_and_genes$'RH-15-Me') / (counts_and_genes$'RH-15-Me' + counts_and_genes$'RH-15-Un')
counts_and_genes$methylation_proportion_16 <- (counts_and_genes$'RH-16-Me') / (counts_and_genes$'RH-16-Me' + counts_and_genes$'RH-16-Un')
counts_and_genes$methylation_proportion_18 <- (counts_and_genes$'RH-18-Me') / (counts_and_genes$'RH-18-Me' + counts_and_genes$'RH-18-Un')
counts_and_genes$methylation_proportion_19 <- (counts_and_genes$'RH-19-Me') / (counts_and_genes$'RH-19-Me' + counts_and_genes$'RH-19-Un')
counts_and_genes$methylation_proportion_21 <- (counts_and_genes$'RH-21-Me') / (counts_and_genes$'RH-21-Me' + counts_and_genes$'RH-21-Un')
counts_and_genes$methylation_proportion_22 <- (counts_and_genes$'RH-22-Me') / (counts_and_genes$'RH-22-Me' + counts_and_genes$'RH-22-Un')

# Set up design matrix, estimate dispersion and run glmFit; then run LRT for contrasts of intrest 

# Group var has your 6 groups - we use a one-way layout for simplicity and given our research questions 
targets_2$group <- gsub(" ", "", targets_2$group, fixed = TRUE)
#targets_2$group<-relevel(as.factor(targets_2$group), ref='SMV normal')
design <- modelMatrixMeth(model.matrix(~ 0 + group, data=targets_2))
design
# targets_2$diet <- relevel(as.factor(targets_2$diet), ref = 'Normal diet')
# targets_2$tissue <- relevel(as.factor(targets_2$tissue), ref = 'Normal myocardium')
# targets_2$surgery <- relevel(as.factor(targets_2$surgery), ref = 'Sham surgery')
# design <- modelMatrixMeth(model.matrix(~ diet + tissue + surgery + diet:tissue + diet:surgery + tissue:surgery, data=targets_2))
#design <- modelMatrixMeth(model.matrix(~0 + group, data = targets_2))

# Comments on dispersion:
  
  "We estimate the NB dispersion for each CpG site using the estimateDisp function. The mean-dispersion relationship of BS-seq data has been studied in the past and no apparent mean-dispersion trend was observed [10]. Therefore, we would not consider a mean-dependent dispersion trend for BS-seq methylation data."

# This is why they use: estimateDisp(trend="none") as dispersion for glmLRT.

# Estimate dispersion and model 
y_glm_disp <- estimateDisp(y, design=design,  trend="none")
y_glm_fit <- glmFit(y_glm_disp, design)

# Contrasts of interest for our study 

# Research Question 1 - SMV normal vs SMV ischemic 
SMVisch_vs_SMVnormal <- glmLRT(y_glm_fit, contrast = makeContrasts(groupSMVnormal - groupSMVischemic, levels=design)) 
SMVisch_vs_SMVnormal

# Research Question 2 - HSMV normal vs HSMV ischemic 
HSMVisch_vs_HSMVnormal <- glmLRT(y_glm_fit, contrast = makeContrasts(groupHSMVnormal - groupHSMVischemic, levels=design)) 
HSMVisch_vs_HSMVnormal

# Research Question 3 - SMV normal vs SMV ischemic vs HSMV normal vs HSMV ischemic 
diffisch_vs_diffnormal <- glmLRT(y_glm_fit, contrast = makeContrasts( (groupSMVnormal - groupSMVischemic) - (groupHSMVnormal - groupHSMVischemic), levels=design)) 
diffisch_vs_diffnormal

# Research Question 4 - MVM vs SMV ischemic 
MVM_vs_SMVischemic <- glmLRT(y_glm_fit, contrast = makeContrasts(groupMVM - groupSMVischemic, levels=design)) 
MVM_vs_SMVischemic

# Research Question 5 - HVM vs HSMV ischemic 
HVM_vs_HSMVischemic <- glmLRT(y_glm_fit, contrast = makeContrasts(groupHVM - groupHSMVischemic, levels=design)) 
HVM_vs_HSMVischemic

# Add corrected p-value columns

# Research Question 1 - SMV normal vs SMV ischemic 
SMVisch_vs_SMVnormal$table$BH <- p.adjust(SMVisch_vs_SMVnormal$table$PValue, method = "BH")
SMVisch_vs_SMVnormal$table$bonferroni <- p.adjust(SMVisch_vs_SMVnormal$table$PValue, method = "bonferroni")
SMVisch_vs_SMVnormal

# Research Question 2 - HSMV normal vs HSMV ischemic 
HSMVisch_vs_HSMVnormal$table$BH <- p.adjust(HSMVisch_vs_HSMVnormal$table$PValue, method = "BH")
HSMVisch_vs_HSMVnormal$table$bonferroni <- p.adjust(HSMVisch_vs_HSMVnormal$table$PValue, method = "bonferroni")
HSMVisch_vs_HSMVnormal 

# Research Question 3 - SMV normal vs SMV ischemic vs HSMV normal vs HSMV ischemic 
diffisch_vs_diffnormal$table$BH <- p.adjust(diffisch_vs_diffnormal$table$PValue, method = "BH")
diffisch_vs_diffnormal$table$bonferroni <- p.adjust(diffisch_vs_diffnormal$table$PValue, method = "bonferroni")
diffisch_vs_diffnormal

# Research Question 4 - MVM vs SMV ischemic 
MVM_vs_SMVischemic$table$BH <- p.adjust(MVM_vs_SMVischemic$table$PValue, method = "BH")
MVM_vs_SMVischemic$table$bonferroni <- p.adjust(MVM_vs_SMVischemic$table$PValue, method = "bonferroni")
MVM_vs_SMVischemic 

# Research Question 5 - HVM vs HSMV ischemic 
HVM_vs_HSMVischemic$table$BH <- p.adjust(HVM_vs_HSMVischemic$table$PValue, method = "BH")
HVM_vs_HSMVischemic$table$bonferroni <- p.adjust(HVM_vs_HSMVischemic$table$PValue, method = "bonferroni")
HVM_vs_HSMVischemic

#################
# Make tables 
#################

# Question 1 
SMVisch_vs_SMVnormal_results <- data.frame(SMVisch_vs_SMVnormal$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
SMVisch_vs_SMVnormal_results <- left_join(SMVisch_vs_SMVnormal_results, counts_and_genes, by = 'loci')
# filter results based on on nominal p-value
SMVisch_vs_SMVnormal_results_filtered <- SMVisch_vs_SMVnormal_results %>% dplyr::filter(PValue < 0.01)
nrow(SMVisch_vs_SMVnormal_results_filtered)
SMVisch_vs_SMVnormal_results_filtered

# Question 2 
HSMVisch_vs_HSMVnormal_results <- data.frame(HSMVisch_vs_HSMVnormal$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
HSMVisch_vs_HSMVnormal_results <- left_join(HSMVisch_vs_HSMVnormal_results, counts_and_genes, by = 'loci')
# filter results based on on nominal p-value
HSMVisch_vs_HSMVnormal_results_filtered <- HSMVisch_vs_HSMVnormal_results %>% dplyr::filter(PValue < 0.01)
nrow(HSMVisch_vs_HSMVnormal_results_filtered)
HSMVisch_vs_HSMVnormal_results_filtered

# Question 3 
diffisch_vs_diffnormal_results <- data.frame(diffisch_vs_diffnormal$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
diffisch_vs_diffnormal_results <- left_join(diffisch_vs_diffnormal_results, counts_and_genes, by = 'loci')
# filter results based on on nominal p-value
diffisch_vs_diffnormal_results_filtered <- diffisch_vs_diffnormal_results %>% dplyr::filter(PValue < 0.01)
nrow(diffisch_vs_diffnormal_results_filtered)
diffisch_vs_diffnormal_results_filtered

# Question 4 
MVM_vs_SMVischemic_results <- data.frame(MVM_vs_SMVischemic$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
MVM_vs_SMVischemic_results <- left_join(MVM_vs_SMVischemic_results, counts_and_genes, by = 'loci')
# filter results based on on nominal p-value
MVM_vs_SMVischemic_results_filtered <- MVM_vs_SMVischemic_results %>% dplyr::filter(PValue < 0.01)
nrow(MVM_vs_SMVischemic_results_filtered)
MVM_vs_SMVischemic_results_filtered

# Question 5 
HVM_vs_HSMVischemic_results <- data.frame(HVM_vs_HSMVischemic$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
HVM_vs_HSMVischemic_results <- left_join(HVM_vs_HSMVischemic_results, counts_and_genes, by = 'loci')
# filter results based on on nominal p-value
HVM_vs_HSMVischemic_results_filtered <- HVM_vs_HSMVischemic_results %>% dplyr::filter(PValue < 0.01)
nrow(HVM_vs_HSMVischemic_results_filtered)
HVM_vs_HSMVischemic_results_filtered

# Print results 
# Full
SMVisch_vs_SMVnormal_results
HSMVisch_vs_HSMVnormal_results
diffisch_vs_diffnormal_results
MVM_vs_SMVischemic_results
HVM_vs_HSMVischemic_results
# Filtered
SMVisch_vs_SMVnormal_results_filtered
HSMVisch_vs_HSMVnormal_results_filtered
diffisch_vs_diffnormal_results_filtered
MVM_vs_SMVischemic_results_filtered
HVM_vs_HSMVischemic_results_filtered


# Write these to excel file 



















###############################################################
#' #write a test output sheet to check logfc, etc.
#' write.table(nothing_across_estrogen_results_additive, '~/nothing_across_estrogen_results_additive.txt', row.names = F)
#' 
#' #make a list of sheets..
#' sheet_list <- list(
#'   #'Est in negative control' = nothing_across_estrogen_results_additive,
#'   #'Est in talc' = talc_across_estrogen_results_additive,
#'   #'Est in titanium' = titanium_across_estrogen_results_additive,
#'   'Talc vs control in Est' = talc_within_estrogen_results_additive,
#'   'Titanium vs control in Est' = titanium_within_estrogen_results_additive,
#'   'Talc vs control in NoEst' = talc_within_noestrogen_results_additive,
#'   'Titanium vs control in NoEst' = titanium_within_noestrogen_results_additive,
#'   'Talc in Est vs Talc in NoEst' = talc_within_estrogen_VS_talc_within_no_estrogen_results_interaction,
#'   'Titanium in Est vs Titanium in NoEst' = titanium_within_estrogen_VS_titanium_within_no_estrogen_results_interaction)
#' 
#' openxlsx::write.xlsx(sheet_list, file = 'DML_stats.xlsx')
#' 
#' nrow(talc_within_estrogen_results_additive)
#' nrow(talc_within_noestrogen_results_additive)
#' nrow(talc_within_estrogen_VS_talc_within_no_estrogen_results_interaction)
#' 
#' nrow(inner_join(talc_within_estrogen_results_additive, talc_within_noestrogen_results_additive, by = 'loci'))
#' nrow(inner_join(talc_within_estrogen_results_additive, talc_within_estrogen_VS_talc_within_no_estrogen_results_interaction, by = 'loci'))
#' nrow(inner_join(talc_within_noestrogen_results_additive, talc_within_estrogen_VS_talc_within_no_estrogen_results_interaction, by = 'loci'))
#' 
#' head(semi_join(talc_within_estrogen_results_additive, talc_within_noestrogen_results_additive, by = 'loci'))
#' 
#' #import the RNAseq file:
#' rnaseq_data <- read.table('filt20_annot_LIMMA_TEvsVeh_005_Metacoreinput_AVG.txt', sep = '\t', header = T)


