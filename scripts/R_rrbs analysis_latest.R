#######################
# title: "RRBS"
# author: "Jordan"
# date: "8/2/2021"
# output: html_document
#######################


# Notes: Do liberal, moderate, and conservative approaches - 3 analyses with results 

# Look back at commands used to download pig genome and look at organism and version and use this with GenomicFeatures package to make txdb - are the genes I care about diff methylated 

# Once you have differentially methylated loci, we will figure out how to annotate them 

# Add exploratory data analysis using MDS 

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
targets_2 <- read_excel("targets.xlsx")
names(targets_2)
targets_2<-as.data.frame(targets_2[,-c(3,4,8)])
head(targets_2)
# Reference target file 
# targets_jos <- read.table('targets_2.txt', sep = '\t', header = T) #targets and targets_2 change additive and control factors

# import those files using readBismark2DGE (can be .cov or .cov.gz)
files <- targets_2$file #these are the names of the cov files 
setwd("/Users/jordan/Desktop/rrbs_r/cov_files")
yall <- readBismark2DGE(files, sample.names=targets_2$sample)

yall <- yall[yall$genes[,'Chr'] == 1 | yall$genes[,'Chr'] == 2| yall$genes[,'Chr'] == 3| yall$genes[,'Chr'] == 4| yall$genes[,'Chr'] == 5| yall$genes[,'Chr'] == 6| yall$genes[,'Chr'] == 7| yall$genes[,'Chr'] == 8| yall$genes[,'Chr'] == 9| yall$genes[,'Chr'] == 10| yall$genes[,'Chr'] == 11| yall$genes[,'Chr'] == 12| yall$genes[,'Chr'] == 13| yall$genes[,'Chr'] == 14| yall$genes[,'Chr'] == 15| yall$genes[,'Chr'] == 16| yall$genes[,'Chr'] == 17| yall$genes[,'Chr'] == 18| yall$genes[,'Chr'] == 'X'| yall$genes[,'Chr'] == 'Y',]

sort_order <- c(1:18, 'X', 'Y')

yall$genes$Chr <- factor(yall$genes$Chr, levels=sort_order)

o <- order(yall$genes$Chr, yall$genes$Locus)

yall <- yall[o,]        

table(yall$genes$Chr)

# Gene filtering and normalization 

# make factor levels of Me and Un repeated along the length of the columns 
Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
Me <- yall$counts[, Methylation=="Me"]
Un <- yall$counts[, Methylation=="Un"]

# combine methylated and unmethylated reads to get total coverage
Coverage <- Me + Un

# Get a list of CpGs with at least 5 counts (methylated and unmethylated) in 3 samples 

HasCoverage <- rowSums(Coverage >= 5) == 3 # requires that sum across methylated and unmethylated counts is at least 60 across three samples 

# List of CpGs that are not 0 both
# this filters out CpG that are never methylated or always methylated as this provides no info about differential methylation 
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0 

# Check sample size after filtering requirements are met 
table(HasCoverage, HasBoth)

# Check sample size after filtering requirements are met 
table(HasCoverage, HasBoth)

# apply those coverage filters
y <- yall[HasCoverage & HasBoth,, keep.lib.sizes=FALSE] 

# assign library sizes to each sample to be a sum of methylated and unmethylated counts
TotalLibSize <- y$samples$lib.size[Methylation=="Me"] + y$samples$lib.size[Methylation=="Un"]
y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples

table(y$genes$Chr)

# make some output tables for later..
library(dplyr)

counts_and_genes <- inner_join(tibble::rownames_to_column(as.data.frame(y$counts), var = 'loci'), tibble::rownames_to_column(as.data.frame(y$genes), var = 'loci'))

counts_and_genes$methylation_proportion_1 <- (counts_and_genes$'1-Me') / (counts_and_genes$'1-Me' + counts_and_genes$'1-Un')
counts_and_genes$methylation_proportion_2 <- (counts_and_genes$'2-Me') / (counts_and_genes$'2-Me' + counts_and_genes$'2-Un')
counts_and_genes$methylation_proportion_3 <- (counts_and_genes$'3-Me') / (counts_and_genes$'3-Me' + counts_and_genes$'3-Un')
counts_and_genes$methylation_proportion_4 <- (counts_and_genes$'4-Me') / (counts_and_genes$'4-Me' + counts_and_genes$'4-Un')
counts_and_genes$methylation_proportion_5 <- (counts_and_genes$'5-Me') / (counts_and_genes$'5-Me' + counts_and_genes$'5-Un')
counts_and_genes$methylation_proportion_6 <- (counts_and_genes$'6-Me') / (counts_and_genes$'6-Me' + counts_and_genes$'6-Un')
counts_and_genes$methylation_proportion_7 <- (counts_and_genes$'7-Me') / (counts_and_genes$'7-Me' + counts_and_genes$'7-Un')
counts_and_genes$methylation_proportion_8 <- (counts_and_genes$'8-Me') / (counts_and_genes$'8-Me' + counts_and_genes$'8-Un')
counts_and_genes$methylation_proportion_10 <- (counts_and_genes$'10-Me') / (counts_and_genes$'10-Me' + counts_and_genes$'10-Un')
counts_and_genes$methylation_proportion_11 <- (counts_and_genes$'11-Me') / (counts_and_genes$'11-Me' + counts_and_genes$'11-Un')
counts_and_genes$methylation_proportion_12 <- (counts_and_genes$'12-Me') / (counts_and_genes$'12-Me' + counts_and_genes$'12-Un')
counts_and_genes$methylation_proportion_14 <- (counts_and_genes$'14-Me') / (counts_and_genes$'14-Me' + counts_and_genes$'14-Un')
counts_and_genes$methylation_proportion_15 <- (counts_and_genes$'15-Me') / (counts_and_genes$'15-Me' + counts_and_genes$'15-Un')
counts_and_genes$methylation_proportion_16 <- (counts_and_genes$'16-Me') / (counts_and_genes$'16-Me' + counts_and_genes$'16-Un')
counts_and_genes$methylation_proportion_18 <- (counts_and_genes$'18-Me') / (counts_and_genes$'18-Me' + counts_and_genes$'18-Un')
counts_and_genes$methylation_proportion_19 <- (counts_and_genes$'19-Me') / (counts_and_genes$'19-Me' + counts_and_genes$'19-Un')
counts_and_genes$methylation_proportion_21 <- (counts_and_genes$'21-Me') / (counts_and_genes$'21-Me' + counts_and_genes$'21-Un')
counts_and_genes$methylation_proportion_22 <- (counts_and_genes$'22-Me') / (counts_and_genes$'22-Me' + counts_and_genes$'22-Un')

# Group var has your 6 groups - we use a one-way layout for simplicity and given our research questions 

targets$groups <- c(
  'groupSMVnormal','groupSMVnormal', 'groupSMVnormal', 
  'groupSMVischemic','groupSMVischemic','groupSMVischemic',
  'MVM', 'MVM', 'MVM',
  'HVM', 'HVM', 'HVM',
  'groupHSMVnormal', 'groupHSMVnormal', 'groupHSMVnormal',
  'groupHSMVischemic', 'groupHSMVischemic', 'groupHSMVischemic')

targets$groups <- gsub(" ", "", targets$groups, fixed = TRUE)

design <- modelMatrixMeth(model.matrix(~ 0 + groups, data=targets))
design

# Estimate dispersion and model 
y_glm_disp <- estimateDisp(y, design=design,  trend="none")
y_glm_fit <- glmFit(y_glm_disp, design)

# Contrasts of interest for our study 

# Research Question 1 - SMV normal vs SMV ischemic 
SMVisch_vs_SMVnormal <- glmLRT(y_glm_fit, contrast = makeContrasts(groupsgroupSMVnormal - groupsgroupSMVischemic, levels=design)) 
SMVisch_vs_SMVnormal

# Research Question 2 - HSMV normal vs HSMV ischemic 
HSMVisch_vs_HSMVnormal <- glmLRT(y_glm_fit, contrast = makeContrasts(groupsgroupHSMVnormal - groupsgroupHSMVischemic, levels=design)) 
HSMVisch_vs_HSMVnormal

# Research Question 3 - SMV normal vs SMV ischemic vs HSMV normal vs HSMV ischemic 
#diffisch_vs_diffnormal <- glmLRT(y_glm_fit, contrast = makeContrasts( (groupsgroupSMVnormal - groupsgroupSMVischemic) - (groupsgroupHSMVnormal - groupsgroupHSMVischemic), levels=design)) 
#diffisch_vs_diffnormal

# Research Question 4 - MVM vs SMV ischemic 
MVM_vs_SMVischemic <- glmLRT(y_glm_fit, contrast = makeContrasts(groupsMVM- groupsgroupSMVischemic, levels=design)) 
MVM_vs_SMVischemic

# Research Question 5 - HVM vs HSMV ischemic 
HVM_vs_HSMVischemic <- glmLRT(y_glm_fit, contrast = makeContrasts(groupsHVM - groupsgroupHSMVischemic, levels=design)) 
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

#Get some annotations

library("biomaRt")
ss_mart <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")

attributes <- listAttributes(ss_mart)
attributes[grep("symbol", attributes$description, ignore.case = TRUE), ]
attributes[grep("transcription", attributes$description, ignore.case = TRUE), ]

ss_tss <- getBM(attributes = c("transcription_start_site", 
                               "chromosome_name", 
                               "start_position",
                               "end_position",
                               "ensembl_gene_id", 
                               "strand", 
                               "uniprot_gn_symbol"), 
                mart = ss_mart)

ss_tss <- as.data.frame(ss_tss)
ss_tss$transcription_end_site <- (ss_tss$transcription_start_site + 1)
ss_tss_gr <- makeGRangesFromDataFrame(ss_tss, start.field = 'transcription_start_site', end.field = 'transcription_end_site', ignore.strand = T)

ss_tss_gr$uniprot_gn_symbol <- ss_tss$uniprot_gn_symbol
ss_tss_gr$ensembl_gene_id <- ss_tss$ensembl_gene_id
ss_tss_gr$gene_start <- ss_tss$start_position
ss_tss_gr$gene_end <- ss_tss$end_position

#get all the genes used in  the DML analysis..
y_loc <- as.data.frame(y$genes)
y_loc$start <- y_loc$Locus
y_loc$end <- (y_loc$start + 1)
row.names(y_loc) <- NULL

y_gr <- makeGRangesFromDataFrame(y_loc, 
                                 ignore.strand = T)

#y_loc and y_gr are the methylation sites
#ss_tss_gr is all the TSS and associated gene symbols, ids, and ranges

#use nearest function in GenomicRanges
y_nearest_tss <- nearest(y = y_gr,
                         x = ss_tss_gr)

dtn <- as.data.frame(distanceToNearest(x = y_gr, subject = ss_tss_gr, ignore.strand=T))

y_gr$distance <- dtn$distance
y_gr$subjectHits <- dtn$subjectHits

#the distance function returns the row indices, use the slice function to pull them out of the y_gr object
tss_slice <- as.data.frame(ss_tss_gr) %>% slice(y_gr$subjectHits)
y_gr_df <- as.data.frame(y_gr)

#add other info like gene symbol, etc.
y_gr_df$TSS_start <- tss_slice$start
y_gr_df$distance <- tss_slice$distance
y_gr_df$uniprot_gn_symbol <- tss_slice$uniprot_gn_symbol
y_gr_df$ensembl_gene_id <- tss_slice$ensembl_gene_id
y_gr_df$gene_start <- tss_slice$gene_start
y_gr_df$gene_end <- tss_slice$gene_end

#make a new 'locus' column to match whats in the results tables
y_gr_df <- unite(y_gr_df, col = 'loci', sep = '-', seqnames:start)

rownames(y_gr_df) <- y_gr_df$locus


#################
# Make tables 
#################

y_gr_df$end <- NULL
y_gr_df$width <- NULL
y_gr_df$strand <- NULL
y_gr_df$subjectHits <- NULL


# Question 1 
SMVisch_vs_SMVnormal_results <- data.frame(SMVisch_vs_SMVnormal$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
SMVisch_vs_SMVnormal_results <- left_join(SMVisch_vs_SMVnormal_results, counts_and_genes, by = 'loci')
#add annotations from `y_gr_df`
SMVisch_vs_SMVnormal_results <- left_join(SMVisch_vs_SMVnormal_results, y_gr_df, by = c('loci' = 'locus'))
# filter results based on on nominal p-value
SMVisch_vs_SMVnormal_results_filtered <- SMVisch_vs_SMVnormal_results %>% dplyr::filter(PValue < 0.01)
nrow(SMVisch_vs_SMVnormal_results_filtered)
SMVisch_vs_SMVnormal_results_filtered

# Question 2 
HSMVisch_vs_HSMVnormal_results <- data.frame(HSMVisch_vs_HSMVnormal$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
HSMVisch_vs_HSMVnormal_results <- left_join(HSMVisch_vs_HSMVnormal_results, counts_and_genes, by = 'loci')
#add annotations from `y_gr_df`
HSMVisch_vs_HSMVnormal_results <- left_join(HSMVisch_vs_HSMVnormal_results, y_gr_df, by = c('loci' = 'locus'))
# filter results based on on nominal p-value
HSMVisch_vs_HSMVnormal_results_filtered <- HSMVisch_vs_HSMVnormal_results %>% dplyr::filter(PValue < 0.01)
nrow(HSMVisch_vs_HSMVnormal_results_filtered)
HSMVisch_vs_HSMVnormal_results_filtered

# # Question 3 
# diffisch_vs_diffnormal_results <- data.frame(diffisch_vs_diffnormal$table) %>% tibble::rownames_to_column(var = 'loci')
# # add counts
# diffisch_vs_diffnormal_results <- left_join(diffisch_vs_diffnormal_results, counts_and_genes, by = 'loci')
# #add annotations from `y_gr_df`
# diffisch_vs_diffnormal_results <- left_join(diffisch_vs_diffnormal_results, y_gr_df, by = c('loci' = 'locus'))
# # filter results based on on nominal p-value
# diffisch_vs_diffnormal_results_filtered <- diffisch_vs_diffnormal_results %>% dplyr::filter(PValue < 0.01)
# nrow(diffisch_vs_diffnormal_results_filtered)
# diffisch_vs_diffnormal_results_filtered

# Question 4 
MVM_vs_SMVischemic_results <- data.frame(MVM_vs_SMVischemic$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
MVM_vs_SMVischemic_results <- left_join(MVM_vs_SMVischemic_results, counts_and_genes, by = 'loci')
#add annotations from `y_gr_df`
MVM_vs_SMVischemic_results <- left_join(MVM_vs_SMVischemic_results, y_gr_df, by = c('loci' = 'locus'))
# filter results based on on nominal p-value
MVM_vs_SMVischemic_results_filtered <- MVM_vs_SMVischemic_results %>% dplyr::filter(PValue < 0.01)
nrow(MVM_vs_SMVischemic_results_filtered)
MVM_vs_SMVischemic_results_filtered

# Question 5 
HVM_vs_HSMVischemic_results <- data.frame(HVM_vs_HSMVischemic$table) %>% tibble::rownames_to_column(var = 'loci')
# add counts
HVM_vs_HSMVischemic_results <- left_join(HVM_vs_HSMVischemic_results, counts_and_genes, by = 'loci')
#add annotations from `y_gr_df`
HVM_vs_HSMVischemic_results <- left_join(HVM_vs_HSMVischemic_results, y_gr_df, by = c('loci' = 'locus'))
# filter results based on on nominal p-value
HVM_vs_HSMVischemic_results_filtered <- HVM_vs_HSMVischemic_results %>% dplyr::filter(PValue < 0.01)
nrow(HVM_vs_HSMVischemic_results_filtered)
HVM_vs_HSMVischemic_results_filtered



# Write worksheets/excel


#First, make the readme tab:
readme_glm <- "glmFit, glmLRT, and makeContrasts functions in edgeR were used to fit a genewise negative binomial GLM and make comparisons of interest based on the research question."
readme_loci <- "`loci` column is the location of a given genomic location. The first number or letter is the chromosome, followed by a `-` symbol and a second number indicates the position or location on the chromosome."
readme_logFC <- "`logFC` indicates the log2-fold change of expression between the two conditions tested. For `SMVisch_vs_SMVnormal_results_filtered` results, positive values indicate higher methylation rates in SMVnormal relative to SMVisch, negative values indicate lower methylation rates in SMVnormal vs SMVisch. For `HSMVisch_vs_HSMVnormal_results_filtered` results, positive values indicate higher methylation rates in HSMVnormal relative to HSMVisch, negative values indicate lower methylation rates in HSMVnormal relative to HSMVisch. For `diffisch_vs_diffnormal_results_filtered` results, positive values indicate ...? For `MVM_vs_SMVischemic_results_filtered` results, positive values indicate higher methylation rates in MVM relative to SMVischemic, negative values indicate lower methylation rates in MVM relative to SMVischemic."
readme_logCPM <- "`logCPM` is the average log2-counts per million, the average taken over all libraries used in the experiment."
readme_LR <- "`LR` is the likelihood ratio statistics (larger LR, smaller p-value"
readme_PValue <- "`PValue` probability of obtaining results at least as extreme as the results actually observed, under the assumption that the null hypothesis is correct"
readme_BH <- "`BH` is Benjamini & Hochberg corrected p-values"
readme_bonferroni <- "`bonferroni` Bonferroni corrected p-values"
readme_MeUn <- "The columns with the format `#-Me` or `#-Un` indicate the methylated and unmethylated counts at a given locus for a particular sample."
readme_Chr <- "`Chr` is the chromosome of the locus"
readme_Locus <- "`Locus` is the position of the locus on the chromosome indicated in the `Chr` column."
readme_methylation_proportion <- "The columns with the format `methylated_proportion_#` indicate the proportion of methylated counts at a given loci for a particular sample, calculated as `methylated counts / (unmethylated counts + methylated counts)` for each sample and loci."
readme_TSS_start <- "`TSS_start` indicates the TSS location nearest to the methylation loci shown."
readme_uniprot_gn_symbol <- "`uniprot_gn_symbol` is the gene symbol (if available) associated with each TSS"
readme_ensembl_gene_id <- "`ensembl_gene_id` is the Ensembl gene ID associated with each TSS"
readme_ensembl_gene_start <- "`gene_start` and `gene_end` are the start and end genomic positions for each gene."

readme_sheet <- rbind(readme_glm, readme_loci, readme_logFC, readme_logCPM, readme_LR, readme_PValue, readme_BH, readme_bonferroni, readme_MeUn, readme_Chr, readme_Locus, readme_methylation_proportion, readme_TSS_start, readme_uniprot_gn_symbol, readme_ensembl_gene_id, readme_ensembl_gene_start)

list_1 <- list(readme_sheet, SMVisch_vs_SMVnormal_results_filtered)
names(list_1) <- c("readme", "SMVisch_vs_SMVnormal")
write.xlsx(list_1, file = '/home/rstudio/dat/SMVisch_vs_SMVnormal_results_filtered.xlsx')

list_2 <- list(readme_sheet, HSMVisch_vs_HSMVnormal_results_filtered)
names(list_2) <- c("readme", "HSMVisch_vs_HSMVnormal")
write.xlsx(list_2, file = '/home/rstudio/dat/HSMVisch_vs_HSMVnormal_results_filtered.xlsx')

# list_3 <- list(readme_sheet, diffisch_vs_diffnormal_results_filtered)
# names(list_3) <- c("readme", "diffisch_vs_diffnormal")
# write.xlsx(list_3, file = '/home/rstudio/dat/diffisch_vs_diffnormal_results_filtered.xlsx')

list_4 <- list(readme_sheet, HVM_vs_HSMVischemic_results_filtered)
names(list_4) <- c("readme", "HVM_vs_HSMVischemic")
write.xlsx(list_4, file = '/home/rstudio/dat/HVM_vs_HSMVischemic_results_filtered.xlsx')

