#### set up ####

setwd(dir = "D:/siobh/Documents/Uni/PhD/Y1/Clozapine Project/Git/clozapine-project/HLA/")
set.seed(666)

#BiocManager::install("HIBAG")

library(HIBAG)
library(parallel)
library(stringr)
library(lmerTest)
library(lme4)
library(tidyverse)
library(datawizard)
library(fastDummies)
library(proxy)


############# genotype imputation with HIBAG ###########
impute <- FALSE

if (impute == TRUE) {
  print("Imputing HLA genotypes using HIBAG")
  # Load the published parameter estimates from European ancestry
  model.list <- get(load(file = "InfiniumGlobal-European-HLA4-hg19.RData"))
  
  # Import your PLINK BED file
  #
  yourgeno <- hlaBED2Geno(bed.fn="CLOZUK3.gsa.qc4.chr6.bed", fam.fn="CLOZUK3.gsa.qc4.chr6.fam", bim.fn="CLOZUK3.gsa.qc4.chr6.bim")
  summary(yourgeno)
  
  #go through all the HLA classic genes
  classic <- c("A","B", "C", "DPB1", "DQA1", "DQB1", "DRB1")
  antonio <- c("B", "DQB1")
  
  #loop through
  for (i in classic) {
    # HLA imputation 
    hla.id <- i
    model <- hlaModelFromObj(model.list[[hla.id]])
    summary(model)
    # HLA allele frequencies
    cbind(frequency = model$hla.freq)
    # best-guess genotypes and all posterior probabilities
    pred.guess <- predict(model, yourgeno, type="response+prob", match.type="Position")
    HLA_geno <- pred.guess$value
    assign(paste0("HLA_", i, "geno"), HLA_geno)
    HLA_pp <- pred.guess$postprob
    assign(paste0("HLA_", i, "pp"), HLA_pp)
  }
  
  
  beepr::beep(sound = 5)
} else {
  print("HLA genotypes already imputed")
}

#### analyse output - A ####

# load
A <- read.table("HLA_Ageno.txt",header=T,col.names=c("LUIN","A.A1","A.A2","A.prob", "A.matching"))


# Define HLA type (making the dummy variables the new way)
A.A1 <- fastDummies::dummy_cols(A, select_columns = c('A.A1'))
A.A1 <- A.A1[,-2:-5]
colnames(A.A1) <- gsub("A.A1_","",colnames(A.A1))
A.A1 <- column_as_rownames(x = A.A1, var = "LUIN")

A.A2 <- fastDummies::dummy_cols(A, select_columns = c('A.A2'))
A.A2 <- A.A2[,-2:-5]
colnames(A.A2) <- gsub("A.A2_","",colnames(A.A2))
A.A2 <- column_as_rownames(x = A.A2, var = "LUIN")

all.A <- cbind(A.A1, A.A2)

LD.A <- sapply(unique(colnames(all.A)),function(x) rowSums(all.A[, colnames(all.A) == x, drop = FALSE]))

# Calculate allele frequencies
LD.A.MAF <- data.frame(A=colnames(LD.A), MAF=NA)
for( i in LD.A.MAF$A ) {LD.A.MAF$MAF[LD.A.MAF$A == i] <- sum(LD.A[,i])/(nrow(LD.A)*2)}
LD.A.table <- LD.A.MAF
LD.A.table$MAF<- round(LD.A.table$MAF, 4)


# Weight probability (as Levin et al. 2015)
LD.A.pp <- read.table("HLA_App.txt",header=T)

LD.A.w <- LD.A
rownames(LD.A.w) <- A$LUIN
for( i in rownames(LD.A.w) ) {
  for (j in colnames(LD.A.w)) {
    p2 <- LD.A.pp[paste0(j,"/",j),paste0("X",i)]
    p1 <- sum(LD.A.pp[grep(j,rownames(LD.A.pp)),paste0("X",i)])-p2
    LD.A.w[i,j]<- (2*p2)+p1}}

# Restrict to MAF>1%
LD.A.qc <- LD.A[,as.character(LD.A.MAF$A[LD.A.MAF$MAF >= 0.01])]
LD.A.w.qc <- LD.A.w[,as.character(LD.A.MAF$A[LD.A.MAF$MAF >= 0.01])]

# those excluded via QC
excl_A <- setdiff(colnames(LD.A.w), colnames(LD.A.w.qc))


# can remove later
load("../CLOZUK3.fbc.pk.gen.RData")
CLOZUK3.fbc.pk.gen <- CLOZUK3.fbc.pk.gen %>%
  mutate(std_age2 = as.vector(scale(age_at_assay,center = T, scale = F))^2) 
CLOZUK3.fbc.pk.std <- datawizard::standardise(x = CLOZUK3.fbc.pk.gen, select = c("dailydose", "age_at_assay", "std_age2", "TDS", "clozapine", "norclozapine"))

# Add covariates included in original model
LD.COV <- CLOZUK3.fbc.pk.std %>% dplyr::select(c("LUIN", "neut_num", "dailydose", "clozapine", "norclozapine", "SEX", "age_at_assay", "std_age2", "TDS"))
CLOZ_IDS <- A$LUIN 
LD.COV <- subset(LD.COV, (LUIN %in% CLOZ_IDS)) # we only want data on the individuals for whom we have our hla data
LD.COV$LUIN <- as.character(LD.COV$LUIN)
LD.Adf <- as.data.frame(LD.A) %>% rownames_as_column(var = 'LUIN') # create column from luin rownames so can merge with covariates
LD.A.COV <- left_join(LD.Adf,LD.COV, by = "LUIN") %>% na.omit() # join together so now we have lt data where each indv matched to their hla allele.

# sort .w.qc col
LD.A.w.qc.df <- as.data.frame(LD.A.w.qc) %>% rownames_to_column(var = "LUIN") # we want the pps in the same format (i.e., matched to lt cloz data)
LD.A.w.qc.COV <- as.data.frame(left_join(LD.COV, LD.A.w.qc.df, by = "LUIN") %>% na.omit() %>% dplyr::select(-c(1:9))) # merge so entries replicated when N has multiple clozuk assays/data points, then get rid of clozuk data - this needs to have same n rows as the cov data so it matches up for the regs. 
LD.A.w.qc.std <- as.data.frame(datawizard::standardise(LD.A.w.qc.COV)) # standardise


# Test linear regression models
#summary(lmer(LD.A.COV$neut_num ~ LD.A.COV$dailydose + LD.A.COV$clozapine + LD.A.COV$norclozapine + LD.A.w.qc.COV[,"02:01"] + LD.A.COV$age_at_assay + LD.A.COV$std_age2  + LD.A.COV$SEX + LD.A.COV$TDS + (1|LD.A.COV$LUIN)))

A_lmm <- lapply(1:ncol(LD.A.w.qc.std), function(x) lmer(LD.A.COV$neut_num ~ LD.A.COV$dailydose + LD.A.COV$clozapine + LD.A.COV$norclozapine + LD.A.w.qc.std[,x] + LD.A.COV$age_at_assay + LD.A.COV$std_age2  + LD.A.COV$SEX + LD.A.COV$TDS + (1|LD.A.COV$LUIN)))
A_summaries <- lapply(A_lmm, summary)
#A_summaries_b <- lapply(A_summaries, function(x) x$coefficients[2, c(1,4)]) # None survive Bonferroni

# get the nice info out so dont need to faff with the big summaries list 

A_lmm_brief <- NULL
for (i in 1:ncol(LD.A.w.qc)){
  B <- A_summaries[[i]]$coefficients[5,1]
  SE <- A_summaries[[i]]$coefficients[5,2]  
  P <- A_summaries[[i]]$coefficients[5,5]
  allele <- colnames(LD.A.w.qc.std)[i]
  A_lmm_brief <- rbind(A_lmm_brief, data.frame(Gene = "A", HLA_Allele=allele, Estimate=B, Std_Error = SE, P=P))
}


beepr::beep(sound = 2)


#### analyse output - B ####

# load
B <- read.table("HLA_Bgeno.txt",header=T,col.names=c("LUIN","B.A1","B.A2","B.prob", "B.matching"))


# Define HLA type (making the dummy variables the new way)
B.A1 <- fastDummies::dummy_cols(B, select_columns = c('B.A1'))
B.A1 <- B.A1[,-2:-5]
colnames(B.A1) <- gsub("B.A1_","",colnames(B.A1))
B.A1 <- column_as_rownames(x = B.A1, var = "LUIN")

B.A2 <- fastDummies::dummy_cols(B, select_columns = c('B.A2'))
B.A2 <- B.A2[,-2:-5]
colnames(B.A2) <- gsub("B.A2_","",colnames(B.A2))
B.A2 <- column_as_rownames(x = B.A2, var = "LUIN")

all.B <- cbind(B.A1, B.A2)

LD.B <- sapply(unique(colnames(all.B)),function(x) rowSums(all.B[, colnames(all.B) == x, drop = FALSE]))


# Calculate allele frequencies
LD.B.MAF <- data.frame(B=colnames(LD.B), MAF=NA)
for( i in LD.B.MAF$B ) {LD.B.MAF$MAF[LD.B.MAF$B == i] <- sum(LD.B[,i])/(nrow(LD.B)*2)}
LD.B.table <- LD.B.MAF
LD.B.table$MAF<- round(LD.B.table$MAF, 4)

# Weight probability (as Levin et al. 2015)
LD.B.pp <- read.table("HLA_Bpp.txt",header=T)

LD.B.w <- LD.B
rownames(LD.B.w) <- B$LUIN
for( i in rownames(LD.B.w) ) {
  for (j in colnames(LD.B.w)) {
    p2 <- LD.B.pp[paste0(j,"/",j),paste0("X",i)]
    p1 <- sum(LD.B.pp[grep(j,rownames(LD.B.pp)),paste0("X",i)])-p2
    LD.B.w[i,j]<- (2*p2)+p1}}

# Restrict to MAF>1%
LD.B.qc <- LD.B[,as.character(LD.B.MAF$B[LD.B.MAF$MAF >= 0.01])]
LD.B.w.qc <- LD.B.w[,as.character(LD.B.MAF$B[LD.B.MAF$MAF >= 0.01])]

# those excluded via QC
excl_B <- setdiff(colnames(LD.B.w), colnames(LD.B.w.qc))


# can remove later
load("../CLOZUK3.fbc.pk.RData")
CLOZUK3.fbc.pk <- CLOZUK3.fbc.pk %>%
  mutate(std_age2 = as.vector(scale(age_at_assay,center = T, scale = F))^2) 
CLOZUK3.fbc.pk.std <- datawizard::standardise(x = CLOZUK3.fbc.pk, select = c("dailydose", "age_at_assay", "std_age2", "TDS", "clozapine", "norclozapine"))

# Add covariates included in original model
LD.COV <- CLOZUK3.fbc.pk.std %>% dplyr::select(c("LUIN", "neut_num", "dailydose", "clozapine", "norclozapine", "SEX", "age_at_assay", "std_age2", "TDS"))
CLOZ_IDS <- B$LUIN 
LD.COV <- subset(LD.COV, (LUIN %in% CLOZ_IDS)) # we only want data on the individuals for whom we have our hla data
LD.COV$LUIN <- as.character(LD.COV$LUIN)
LD.Bdf <- as.data.frame(LD.B) %>% rownames_as_column(var = 'LUIN') # create column from luin rownames so can merge with covariates
LD.B.COV <- left_join(LD.Bdf,LD.COV, by = "LUIN") %>% na.omit() # join together so now we have lt data where each indv matched to their hla allele.

# sort .w.qc col
LD.B.w.qc.df <- as.data.frame(LD.B.w.qc) %>% rownames_to_column(var = "LUIN") # we want the pps in the same format (i.e., matched to lt cloz data)
LD.B.w.qc.COV <- as.data.frame(left_join(LD.COV, LD.B.w.qc.df, by = "LUIN") %>% na.omit() %>% dplyr::select(-c(1:9))) # merge so entries replicated when N has multiple clozuk assays/data points, then get rid of clozuk data - this needs to have same n rows as the cov data so it matches up for the regs. 
LD.B.w.qc.std <- as.data.frame(datawizard::standardise(LD.B.w.qc.COV)) # standardise


# Test linear regression models
#summary(lmer(LD.B$neut_num ~ LD.B.COV$dailydose + LD.B.COV$clozapine + LD.B.COV$norclozapine + LD.B.w.qc.COV[,"02:01"] + LD.B.COV$age_at_assay + LD.B.COV$std_age2  + LD.B.COV$SEX + LD.B.COV$TDS + (1|LD.B.COV$LUIN)))

B_lmm <- lapply(1:ncol(LD.B.w.qc.std), function(x) lmer(LD.B.COV$neut_num ~ LD.B.COV$dailydose + LD.B.COV$clozapine + LD.B.COV$norclozapine + LD.B.w.qc.std[,x] + LD.B.COV$age_at_assay + LD.B.COV$std_age2  + LD.B.COV$SEX + LD.B.COV$TDS + (1|LD.B.COV$LUIN)))
B_summaries <- lapply(B_lmm, summary)
#B_summaries_b <- lapply(B_summaries, function(x) x$coefficients[2, c(1,4)]) # None survive Bonferroni

# get the nice info out so dont need to faff with the big summaries list 

B_lmm_brief <- NULL
for (i in 1:ncol(LD.B.w.qc)){
  B <- B_summaries[[i]]$coefficients[5,1]
  SE <- B_summaries[[i]]$coefficients[5,2]  
  P <- B_summaries[[i]]$coefficients[5,5]
  allele <- colnames(LD.B.w.qc.std)[i]
  B_lmm_brief <- rbind(B_lmm_brief, data.frame(Gene = "B", HLA_Allele=allele, Estimate=B, Std_Error = SE, P=P))
}

beepr::beep(sound = 2)




#### analyse output - C ####

# load
C <- read.table("HLA_Cgeno.txt",header=T,col.names=c("LUIN","C.A1","C.A2","C.prob", "C.matching"))


# Define HLA type (making the dummy variables the new way)
C.A1 <- fastDummies::dummy_cols(C, select_columns = c('C.A1'))
C.A1 <- C.A1[,-2:-5]
colnames(C.A1) <- gsub("C.A1_","",colnames(C.A1))
C.A1 <- column_as_rownames(x = C.A1, var = "LUIN")

C.A2 <- fastDummies::dummy_cols(C, select_columns = c('C.A2'))
C.A2 <- C.A2[,-2:-5]
colnames(C.A2) <- gsub("C.A2_","",colnames(C.A2))
C.A2 <- column_as_rownames(x = C.A2, var = "LUIN")

all.C <- cbind(C.A1, C.A2)

LD.C <- sapply(unique(colnames(all.C)),function(x) rowSums(all.C[, colnames(all.C) == x, drop = FALSE]))


# Calculate allele frequencies
LD.C.MAF <- data.frame(C=colnames(LD.C), MAF=NA)
for( i in LD.C.MAF$C ) {LD.C.MAF$MAF[LD.C.MAF$C == i] <- sum(LD.C[,i])/(nrow(LD.C)*2)}
LD.C.table <- LD.C.MAF
LD.C.table$MAF<- round(LD.C.table$MAF, 4)

# Weight probability (as Levin et al. 2015)
LD.C.pp <- read.table("HLA_Cpp.txt",header=T)

LD.C.w <- LD.C
rownames(LD.C.w) <- C$LUIN
for( i in rownames(LD.C.w) ) {
  for (j in colnames(LD.C.w)) {
    p2 <- LD.C.pp[paste0(j,"/",j),paste0("X",i)]
    p1 <- sum(LD.C.pp[grep(j,rownames(LD.C.pp)),paste0("X",i)])-p2
    LD.C.w[i,j]<- (2*p2)+p1}}

# Restrict to MAF>1%
LD.C.qc <- LD.C[,as.character(LD.C.MAF$C[LD.C.MAF$MAF >= 0.01])]
LD.C.w.qc <- LD.C.w[,as.character(LD.C.MAF$C[LD.C.MAF$MAF >= 0.01])]

# those excluded via QC
excl_C <- setdiff(colnames(LD.C.w), colnames(LD.C.w.qc))


# can remove later
load("../CLOZUK3.fbc.pk.RData")
CLOZUK3.fbc.pk <- CLOZUK3.fbc.pk %>%
  mutate(std_age2 = as.vector(scale(age_at_assay,center = T, scale = F))^2) 
CLOZUK3.fbc.pk.std <- datawizard::standardise(x = CLOZUK3.fbc.pk, select = c("dailydose", "age_at_assay", "std_age2", "TDS", "clozapine", "norclozapine"))

# Add covariates included in original model
LD.COV <- CLOZUK3.fbc.pk.std %>% dplyr::select(c("LUIN", "neut_num", "dailydose", "clozapine", "norclozapine", "SEX", "age_at_assay", "std_age2", "TDS"))
CLOZ_IDS <- C$LUIN 
LD.COV <- subset(LD.COV, (LUIN %in% CLOZ_IDS)) # we only want data on the individuals for whom we have our hla data
LD.COV$LUIN <- as.character(LD.COV$LUIN)
LD.Cdf <- as.data.frame(LD.C) %>% rownames_as_column(var = 'LUIN') # create column from luin rownames so can merge with covariates
LD.C.COV <- left_join(LD.Cdf,LD.COV, by = "LUIN") %>% na.omit() # join together so now we have lt data where each indv matched to their hla allele.

# sort .w.qc col
LD.C.w.qc.df <- as.data.frame(LD.C.w.qc) %>% rownames_to_column(var = "LUIN") # we want the pps in the same format (i.e., matched to lt cloz data)
LD.C.w.qc.COV <- as.data.frame(left_join(LD.COV, LD.C.w.qc.df, by = "LUIN") %>% na.omit() %>% dplyr::select(-c(1:9))) # merge so entries replicated when N has multiple clozuk assays/data points, then get rid of clozuk data - this needs to have same n rows as the cov data so it matches up for the regs. 
LD.C.w.qc.std <- as.data.frame(datawizard::standardise(LD.C.w.qc.COV)) # standardise


# Test linear regression models
#summary(lmer(LD.C$neut_num ~ LD.C.COV$dailydose + LD.C.COV$clozapine + LD.C.COV$norclozapine + LD.C.w.qc.COV[,"02:01"] + LD.C.COV$age_at_assay + LD.C.COV$std_age2  + LD.C.COV$SEX + LD.C.COV$TDS + (1|LD.C.COV$LUIN)))

C_lmm <- lapply(1:ncol(LD.C.w.qc.std), function(x) lmer(LD.C.COV$neut_num ~ LD.C.COV$dailydose + LD.C.COV$clozapine + LD.C.COV$norclozapine + LD.C.w.qc.std[,x] + LD.C.COV$age_at_assay + LD.C.COV$std_age2  + LD.C.COV$SEX + LD.C.COV$TDS + (1|LD.C.COV$LUIN)))
C_summaries <- lapply(C_lmm, summary)
#C_summaries_b <- lapply(C_summaries, function(x) x$coefficients[2, c(1,4)]) # None survive Bonferroni

# get the nice info out so dont need to faff with the big summaries list 

C_lmm_brief <- NULL
for (i in 1:ncol(LD.C.w.qc)){
  B <- C_summaries[[i]]$coefficients[5,1]
  SE <- C_summaries[[i]]$coefficients[5,2]  
  P <- C_summaries[[i]]$coefficients[5,5]
  allele <- colnames(LD.C.w.qc.std)[i]
  C_lmm_brief <- rbind(C_lmm_brief, data.frame(Gene = "C", HLA_Allele=allele, Estimate=B, Std_Error = SE, P=P))
}


beepr::beep(sound = 2)


#### analyse output - DPB1 ####

# load
DPB1 <- read.table("HLA_DPB1geno.txt",header=T,col.names=c("LUIN","DPB1.A1","DPB1.A2","DPB1.prob", "DPB1.matching"))


# Define HLA type (making the dummy variables the new way)
DPB1.A1 <- fastDummies::dummy_cols(DPB1, select_columns = c('DPB1.A1'))
DPB1.A1 <- DPB1.A1[,-2:-5]
colnames(DPB1.A1) <- gsub("DPB1.A1_","",colnames(DPB1.A1))
DPB1.A1 <- column_as_rownames(x = DPB1.A1, var = "LUIN")

DPB1.A2 <- fastDummies::dummy_cols(DPB1, select_columns = c('DPB1.A2'))
DPB1.A2 <- DPB1.A2[,-2:-5]
colnames(DPB1.A2) <- gsub("DPB1.A2_","",colnames(DPB1.A2))
DPB1.A2 <- column_as_rownames(x = DPB1.A2, var = "LUIN")

all.DPB1 <- cbind(DPB1.A1, DPB1.A2)

LD.DPB1 <- sapply(unique(colnames(all.DPB1)),function(x) rowSums(all.DPB1[, colnames(all.DPB1) == x, drop = FALSE]))


# Calculate allele frequencies
LD.DPB1.MAF <- data.frame(DPB1=colnames(LD.DPB1), MAF=NA)
for( i in LD.DPB1.MAF$DPB1 ) {LD.DPB1.MAF$MAF[LD.DPB1.MAF$DPB1 == i] <- sum(LD.DPB1[,i])/(nrow(LD.DPB1)*2)}
LD.DPB1.table <- LD.DPB1.MAF
LD.DPB1.table$MAF<- round(LD.DPB1.table$MAF, 4)

# Weight probability (as Levin et al. 2015)
LD.DPB1.pp <- read.table("HLA_DPB1pp.txt",header=T)

LD.DPB1.w <- LD.DPB1
rownames(LD.DPB1.w) <- DPB1$LUIN
for( i in rownames(LD.DPB1.w) ) {
  for (j in colnames(LD.DPB1.w)) {
    p2 <- LD.DPB1.pp[paste0(j,"/",j),paste0("X",i)]
    p1 <- sum(LD.DPB1.pp[grep(j,rownames(LD.DPB1.pp)),paste0("X",i)])-p2
    LD.DPB1.w[i,j]<- (2*p2)+p1}}

# Restrict to MAF>1%
LD.DPB1.qc <- LD.DPB1[,as.character(LD.DPB1.MAF$DPB1[LD.DPB1.MAF$MAF >= 0.01])]
LD.DPB1.w.qc <- LD.DPB1.w[,as.character(LD.DPB1.MAF$DPB1[LD.DPB1.MAF$MAF >= 0.01])]

# those excluded via QC
excl_DPB1 <- setdiff(colnames(LD.DPB1.w), colnames(LD.DPB1.w.qc))


# can remove later
load("../CLOZUK3.fbc.pk.RData")
CLOZUK3.fbc.pk <- CLOZUK3.fbc.pk %>%
  mutate(std_age2 = as.vector(scale(age_at_assay,center = T, scale = F))^2) 
CLOZUK3.fbc.pk.std <- datawizard::standardise(x = CLOZUK3.fbc.pk, select = c("dailydose", "age_at_assay", "std_age2", "TDS", "clozapine", "norclozapine"))

# Add covariates included in original model
LD.COV <- CLOZUK3.fbc.pk.std %>% dplyr::select(c("LUIN", "neut_num", "dailydose", "clozapine", "norclozapine", "SEX", "age_at_assay", "std_age2", "TDS"))
CLOZ_IDS <- DPB1$LUIN 
LD.COV <- subset(LD.COV, (LUIN %in% CLOZ_IDS)) # we only want data on the individuals for whom we have our hla data
LD.COV$LUIN <- as.character(LD.COV$LUIN)
LD.DPB1df <- as.data.frame(LD.DPB1) %>% rownames_as_column(var = 'LUIN') # create column from luin rownames so can merge with covariates
LD.DPB1.COV <- left_join(LD.DPB1df,LD.COV, by = "LUIN") %>% na.omit() # join together so now we have lt data where each indv matched to their hla allele.

# sort .w.qc col
LD.DPB1.w.qc.df <- as.data.frame(LD.DPB1.w.qc) %>% rownames_to_column(var = "LUIN") # we want the pps in the same format (i.e., matched to lt cloz data)
LD.DPB1.w.qc.COV <- as.data.frame(left_join(LD.COV, LD.DPB1.w.qc.df, by = "LUIN") %>% na.omit() %>% dplyr::select(-c(1:9))) # merge so entries replicated when N has multiple clozuk assays/data points, then get rid of clozuk data - this needs to have same n rows as the cov data so it matches up for the regs. 
LD.DPB1.w.qc.std <- as.data.frame(datawizard::standardise(LD.DPB1.w.qc.COV)) # standardise


# Test linear regression models
#summary(lmer(LD.DPB1$neut_num ~ LD.DPB1.COV$dailydose + LD.DPB1.COV$clozapine + LD.DPB1.COV$norclozapine + LD.DPB1.w.qc.COV[,"02:01"] + LD.DPB1.COV$age_at_assay + LD.DPB1.COV$std_age2  + LD.DPB1.COV$SEX + LD.DPB1.COV$TDS + (1|LD.DPB1.COV$LUIN)))

DPB1_lmm <- lapply(1:ncol(LD.DPB1.w.qc.std), function(x) lmer(LD.DPB1.COV$neut_num ~ LD.DPB1.COV$dailydose + LD.DPB1.COV$clozapine + LD.DPB1.COV$norclozapine + LD.DPB1.w.qc.std[,x] + LD.DPB1.COV$age_at_assay + LD.DPB1.COV$std_age2  + LD.DPB1.COV$SEX + LD.DPB1.COV$TDS + (1|LD.DPB1.COV$LUIN)))
DPB1_summaries <- lapply(DPB1_lmm, summary)
#DPB1_summaries_b <- lapply(DPB1_summaries, function(x) x$coefficients[2, c(1,4)]) # None survive Bonferroni

# get the nice info out so dont need to faff with the big summaries list 

DPB1_lmm_brief <- NULL
for (i in 1:ncol(LD.DPB1.w.qc)){
  B <- DPB1_summaries[[i]]$coefficients[5,1]
  SE <- DPB1_summaries[[i]]$coefficients[5,2]  
  P <- DPB1_summaries[[i]]$coefficients[5,5]
  allele <- colnames(LD.DPB1.w.qc.std)[i]
  DPB1_lmm_brief <- rbind(DPB1_lmm_brief, data.frame(Gene = "DPB1", HLA_Allele=allele, Estimate=B, Std_Error = SE, P=P))
}


beepr::beep(sound = 5)




#### analyse output - DQA1 ####

# load
DQA1 <- read.table("HLA_DQA1geno.txt",header=T,col.names=c("LUIN","DQA1.A1","DQA1.A2","DQA1.prob", "DQA1.matching"))


# Define HLA type (making the dummy variables the new way)
DQA1.A1 <- fastDummies::dummy_cols(DQA1, select_columns = c('DQA1.A1'))
DQA1.A1 <- DQA1.A1[,-2:-5]
colnames(DQA1.A1) <- gsub("DQA1.A1_","",colnames(DQA1.A1))
DQA1.A1 <- column_as_rownames(x = DQA1.A1, var = "LUIN")

DQA1.A2 <- fastDummies::dummy_cols(DQA1, select_columns = c('DQA1.A2'))
DQA1.A2 <- DQA1.A2[,-2:-5]
colnames(DQA1.A2) <- gsub("DQA1.A2_","",colnames(DQA1.A2))
DQA1.A2 <- column_as_rownames(x = DQA1.A2, var = "LUIN")

all.DQA1 <- cbind(DQA1.A1, DQA1.A2)

LD.DQA1 <- sapply(unique(colnames(all.DQA1)),function(x) rowSums(all.DQA1[, colnames(all.DQA1) == x, drop = FALSE]))


# Calculate allele frequencies
LD.DQA1.MAF <- data.frame(DQA1=colnames(LD.DQA1), MAF=NA)
for( i in LD.DQA1.MAF$DQA1 ) {LD.DQA1.MAF$MAF[LD.DQA1.MAF$DQA1 == i] <- sum(LD.DQA1[,i])/(nrow(LD.DQA1)*2)}
LD.DQA1.table <- LD.DQA1.MAF
LD.DQA1.table$MAF<- round(LD.DQA1.table$MAF, 4)

# Weight probability (as Levin et al. 2015)
LD.DQA1.pp <- read.table("HLA_DQA1pp.txt",header=T)

LD.DQA1.w <- LD.DQA1
rownames(LD.DQA1.w) <- DQA1$LUIN
for( i in rownames(LD.DQA1.w) ) {
  for (j in colnames(LD.DQA1.w)) {
    p2 <- LD.DQA1.pp[paste0(j,"/",j),paste0("X",i)]
    p1 <- sum(LD.DQA1.pp[grep(j,rownames(LD.DQA1.pp)),paste0("X",i)])-p2
    LD.DQA1.w[i,j]<- (2*p2)+p1}}

# Restrict to MAF>1%
LD.DQA1.qc <- LD.DQA1[,as.character(LD.DQA1.MAF$DQA1[LD.DQA1.MAF$MAF >= 0.01])]
LD.DQA1.w.qc <- LD.DQA1.w[,as.character(LD.DQA1.MAF$DQA1[LD.DQA1.MAF$MAF >= 0.01])]

# those excluded via QC
excl_DQA1 <- setdiff(colnames(LD.DQA1.w), colnames(LD.DQA1.w.qc))


# can remove later
load("../CLOZUK3.fbc.pk.RData")
CLOZUK3.fbc.pk <- CLOZUK3.fbc.pk %>%
  mutate(std_age2 = as.vector(scale(age_at_assay,center = T, scale = F))^2) 
CLOZUK3.fbc.pk.std <- datawizard::standardise(x = CLOZUK3.fbc.pk, select = c("dailydose", "age_at_assay", "std_age2", "TDS", "clozapine", "norclozapine"))

# Add covariates included in original model
LD.COV <- CLOZUK3.fbc.pk.std %>% dplyr::select(c("LUIN", "neut_num", "dailydose", "clozapine", "norclozapine", "SEX", "age_at_assay", "std_age2", "TDS"))
CLOZ_IDS <- DQA1$LUIN 
LD.COV <- subset(LD.COV, (LUIN %in% CLOZ_IDS)) # we only want data on the individuals for whom we have our hla data
LD.COV$LUIN <- as.character(LD.COV$LUIN)
LD.DQA1df <- as.data.frame(LD.DQA1) %>% rownames_as_column(var = 'LUIN') # create column from luin rownames so can merge with covariates
LD.DQA1.COV <- left_join(LD.DQA1df,LD.COV, by = "LUIN") %>% na.omit() # join together so now we have lt data where each indv matched to their hla allele.

# sort .w.qc col
LD.DQA1.w.qc.df <- as.data.frame(LD.DQA1.w.qc) %>% rownames_to_column(var = "LUIN") # we want the pps in the same format (i.e., matched to lt cloz data)
LD.DQA1.w.qc.COV <- as.data.frame(left_join(LD.COV, LD.DQA1.w.qc.df, by = "LUIN") %>% na.omit() %>% dplyr::select(-c(1:9))) # merge so entries replicated when N has multiple clozuk assays/data points, then get rid of clozuk data - this needs to have same n rows as the cov data so it matches up for the regs. 
LD.DQA1.w.qc.std <- as.data.frame(datawizard::standardise(LD.DQA1.w.qc.COV)) # standardise


# Test linear regression models
#summary(lmer(LD.DQA1$neut_num ~ LD.DQA1.COV$dailydose + LD.DQA1.COV$clozapine + LD.DQA1.COV$norclozapine + LD.DQA1.w.qc.COV[,"02:01"] + LD.DQA1.COV$age_at_assay + LD.DQA1.COV$std_age2  + LD.DQA1.COV$SEX + LD.DQA1.COV$TDS + (1|LD.DQA1.COV$LUIN)))

DQA1_lmm <- lapply(1:ncol(LD.DQA1.w.qc.std), function(x) lmer(LD.DQA1.COV$neut_num ~ LD.DQA1.COV$dailydose + LD.DQA1.COV$clozapine + LD.DQA1.COV$norclozapine + LD.DQA1.w.qc.std[,x] + LD.DQA1.COV$age_at_assay + LD.DQA1.COV$std_age2  + LD.DQA1.COV$SEX + LD.DQA1.COV$TDS + (1|LD.DQA1.COV$LUIN)))
DQA1_summaries <- lapply(DQA1_lmm, summary)
#DQA1_summaries_b <- lapply(DQA1_summaries, function(x) x$coefficients[2, c(1,4)]) # None survive Bonferroni

# get the nice info out so dont need to faff with the big summaries list 

DQA1_lmm_brief <- NULL
for (i in 1:ncol(LD.DQA1.w.qc)){
  gene <- DQA1
  B <- DQA1_summaries[[i]]$coefficients[5,1]
  SE <- DQA1_summaries[[i]]$coefficients[5,2]  
  P <- DQA1_summaries[[i]]$coefficients[5,5]
  allele <- colnames(LD.DQA1.w.qc.std)[i]
  DQA1_lmm_brief <- rbind(DQA1_lmm_brief, data.frame(Gene = "DQA1", HLA_Allele=allele, Estimate=B, Std_Error = SE, P=P))
}



beepr::beep(sound = 2)







#### analyse output - DQB1 ####

# load
DQB1 <- read.table("HLA_DQB1geno.txt",header=T,col.names=c("LUIN","DQB1.A1","DQB1.A2","DQB1.prob", "DQB1.matching"))


# Define HLA type (making the dummy variables the new way)
DQB1.A1 <- fastDummies::dummy_cols(DQB1, select_columns = c('DQB1.A1'))
DQB1.A1 <- DQB1.A1[,-2:-5]
colnames(DQB1.A1) <- gsub("DQB1.A1_","",colnames(DQB1.A1))
DQB1.A1 <- column_as_rownames(x = DQB1.A1, var = "LUIN")

DQB1.A2 <- fastDummies::dummy_cols(DQB1, select_columns = c('DQB1.A2'))
DQB1.A2 <- DQB1.A2[,-2:-5]
colnames(DQB1.A2) <- gsub("DQB1.A2_","",colnames(DQB1.A2))
DQB1.A2 <- column_as_rownames(x = DQB1.A2, var = "LUIN")

all.DQB1 <- cbind(DQB1.A1, DQB1.A2)

LD.DQB1 <- sapply(unique(colnames(all.DQB1)),function(x) rowSums(all.DQB1[, colnames(all.DQB1) == x, drop = FALSE]))

# Calculate allele frequencies
LD.DQB1.MAF <- data.frame(DQB1=colnames(LD.DQB1), MAF=NA)
for( i in LD.DQB1.MAF$DQB1 ) {LD.DQB1.MAF$MAF[LD.DQB1.MAF$DQB1 == i] <- sum(LD.DQB1[,i])/(nrow(LD.DQB1)*2)}
LD.DQB1.table <- LD.DQB1.MAF
LD.DQB1.table$MAF<- round(LD.DQB1.table$MAF, 4)


# Weight probability (as Levin et al. 2015)
LD.DQB1.pp <- read.table("HLA_DQB1pp.txt",header=T)

LD.DQB1.w <- LD.DQB1
rownames(LD.DQB1.w) <- DQB1$LUIN
for( i in rownames(LD.DQB1.w) ) {
  for (j in colnames(LD.DQB1.w)) {
    p2 <- LD.DQB1.pp[paste0(j,"/",j),paste0("X",i)]
    p1 <- sum(LD.DQB1.pp[grep(j,rownames(LD.DQB1.pp)),paste0("X",i)])-p2
    LD.DQB1.w[i,j]<- (2*p2)+p1}}

# Restrict to MAF>1%
LD.DQB1.qc <- LD.DQB1[,as.character(LD.DQB1.MAF$DQB1[LD.DQB1.MAF$MAF >= 0.01])]
LD.DQB1.w.qc <- LD.DQB1.w[,as.character(LD.DQB1.MAF$DQB1[LD.DQB1.MAF$MAF >= 0.01])]

# those excluded via QC
excl_DQB1 <- setdiff(colnames(LD.DQB1.w), colnames(LD.DQB1.w.qc))


# can remove later
load("../CLOZUK3.fbc.pk.RData")
CLOZUK3.fbc.pk <- CLOZUK3.fbc.pk %>%
  mutate(std_age2 = as.vector(scale(age_at_assay,center = T, scale = F))^2) 
CLOZUK3.fbc.pk.std <- datawizard::standardise(x = CLOZUK3.fbc.pk, select = c("dailydose", "age_at_assay", "std_age2", "TDS", "clozapine", "norclozapine"))

# Add covariates included in original model
LD.COV <- CLOZUK3.fbc.pk.std %>% dplyr::select(c("LUIN", "neut_num", "dailydose", "clozapine", "norclozapine", "SEX", "age_at_assay", "std_age2", "TDS"))
CLOZ_IDS <- DQB1$LUIN 
LD.COV <- subset(LD.COV, (LUIN %in% CLOZ_IDS)) # we only want data on the individuals for whom we have our hla data
LD.COV$LUIN <- as.character(LD.COV$LUIN)
LD.DQB1df <- as.data.frame(LD.DQB1) %>% rownames_as_column(var = 'LUIN') # create column from luin rownames so can merge with covariates
LD.DQB1.COV <- left_join(LD.DQB1df,LD.COV, by = "LUIN") %>% na.omit() # join together so now we have lt data where each indv matched to their hla allele.

# sort .w.qc col
LD.DQB1.w.qc.df <- as.data.frame(LD.DQB1.w.qc) %>% rownames_to_column(var = "LUIN") # we want the pps in the same format (i.e., matched to lt cloz data)
LD.DQB1.w.qc.COV <- as.data.frame(left_join(LD.COV, LD.DQB1.w.qc.df, by = "LUIN") %>% na.omit() %>% dplyr::select(-c(1:9))) # merge so entries replicated when N has multiple clozuk assays/data points, then get rid of clozuk data - this needs to have same n rows as the cov data so it matches up for the regs. 
LD.DQB1.w.qc.std <- as.data.frame(datawizard::standardise(LD.DQB1.w.qc.COV)) # standardise


# Test linear regression models
#summary(lmer(LD.DQB1$neut_num ~ LD.DQB1.COV$dailydose + LD.DQB1.COV$clozapine + LD.DQB1.COV$norclozapine + LD.DQB1.w.qc.COV[,"02:01"] + LD.DQB1.COV$age_at_assay + LD.DQB1.COV$std_age2  + LD.DQB1.COV$SEX + LD.DQB1.COV$TDS + (1|LD.DQB1.COV$LUIN)))

DQB1_lmm <- lapply(1:ncol(LD.DQB1.w.qc.std), function(x) lmer(LD.DQB1.COV$neut_num ~ LD.DQB1.COV$dailydose + LD.DQB1.COV$clozapine + LD.DQB1.COV$norclozapine + LD.DQB1.w.qc.std[,x] + LD.DQB1.COV$age_at_assay + LD.DQB1.COV$std_age2  + LD.DQB1.COV$SEX + LD.DQB1.COV$TDS + (1|LD.DQB1.COV$LUIN)))
DQB1_summaries <- lapply(DQB1_lmm, summary)
#DQB1_summaries_b <- lapply(DQB1_summaries, function(x) x$coefficients[2, c(1,4)]) # None survive Bonferroni

# get the nice info out so dont need to faff with the big summaries list 

DQB1_lmm_brief <- NULL
for (i in 1:ncol(LD.DQB1.w.qc)){
  B <- DQB1_summaries[[i]]$coefficients[5,1]
  SE <- DQB1_summaries[[i]]$coefficients[5,2]  
  P <- DQB1_summaries[[i]]$coefficients[5,5]
  allele <- colnames(LD.DQB1.w.qc.std)[i]
  DQB1_lmm_brief <- rbind(DQB1_lmm_brief, data.frame(Gene = "DQB1", HLA_Allele=allele, Estimate=B, Std_Error = SE, P=P))
}

# adjust for multiple comparisons which i assume is number of alleles tested for each gene ? or is it number of alleles for all genes?

beepr::beep(sound = 2)







#### analyse output - DRB1 ####

# load
DRB1 <- read.table("HLA_DRB1geno.txt",header=T,col.names=c("LUIN","DRB1.A1","DRB1.A2","DRB1.prob", "DRB1.matching"))


# Define HLA type (making the dummy variables the new way)
DRB1.A1 <- fastDummies::dummy_cols(DRB1, select_columns = c('DRB1.A1'))
DRB1.A1 <- DRB1.A1[,-2:-5]
colnames(DRB1.A1) <- gsub("DRB1.A1_","",colnames(DRB1.A1))
DRB1.A1 <- column_as_rownames(x = DRB1.A1, var = "LUIN")

DRB1.A2 <- fastDummies::dummy_cols(DRB1, select_columns = c('DRB1.A2'))
DRB1.A2 <- DRB1.A2[,-2:-5]
colnames(DRB1.A2) <- gsub("DRB1.A2_","",colnames(DRB1.A2))
DRB1.A2 <- column_as_rownames(x = DRB1.A2, var = "LUIN")

all.DRB1 <- cbind(DRB1.A1, DRB1.A2)

LD.DRB1 <- sapply(unique(colnames(all.DRB1)),function(x) rowSums(all.DRB1[, colnames(all.DRB1) == x, drop = FALSE]))


# Calculate allele frequencies
LD.DRB1.MAF <- data.frame(DRB1=colnames(LD.DRB1), MAF=NA)
for( i in LD.DRB1.MAF$DRB1 ) {LD.DRB1.MAF$MAF[LD.DRB1.MAF$DRB1 == i] <- sum(LD.DRB1[,i])/(nrow(LD.DRB1)*2)}
LD.DRB1.table <- LD.DRB1.MAF
LD.DRB1.table$MAF<- round(LD.DRB1.table$MAF, 4)

# Weight probability (as Levin et al. 2015)
LD.DRB1.pp <- read.table("HLA_DRB1pp.txt",header=T)

LD.DRB1.w <- LD.DRB1
rownames(LD.DRB1.w) <- DRB1$LUIN
for( i in rownames(LD.DRB1.w) ) {
  for (j in colnames(LD.DRB1.w)) {
    p2 <- LD.DRB1.pp[paste0(j,"/",j),paste0("X",i)]
    p1 <- sum(LD.DRB1.pp[grep(j,rownames(LD.DRB1.pp)),paste0("X",i)])-p2
    LD.DRB1.w[i,j]<- (2*p2)+p1}}

# Restrict to MAF>1%
LD.DRB1.qc <- LD.DRB1[,as.character(LD.DRB1.MAF$DRB1[LD.DRB1.MAF$MAF >= 0.01])]
LD.DRB1.w.qc <- LD.DRB1.w[,as.character(LD.DRB1.MAF$DRB1[LD.DRB1.MAF$MAF >= 0.01])]

# those excluded via QC
excl_DRB1 <- setdiff(colnames(LD.DRB1.w), colnames(LD.DRB1.w.qc))


# can remove later
load("../CLOZUK3.fbc.pk.RData")
CLOZUK3.fbc.pk <- CLOZUK3.fbc.pk %>%
  mutate(std_age2 = as.vector(scale(age_at_assay,center = T, scale = F))^2) 
CLOZUK3.fbc.pk.std <- datawizard::standardise(x = CLOZUK3.fbc.pk, select = c("dailydose", "age_at_assay", "std_age2", "TDS", "clozapine", "norclozapine"))

# Add covariates included in original model
LD.COV <- CLOZUK3.fbc.pk.std %>% dplyr::select(c("LUIN", "neut_num", "dailydose", "clozapine", "norclozapine", "SEX", "age_at_assay", "std_age2", "TDS"))
CLOZ_IDS <- DRB1$LUIN 
LD.COV <- subset(LD.COV, (LUIN %in% CLOZ_IDS)) # we only want data on the individuals for whom we have our hla data
LD.COV$LUIN <- as.character(LD.COV$LUIN)
LD.DRB1df <- as.data.frame(LD.DRB1) %>% rownames_as_column(var = 'LUIN') # create column from luin rownames so can merge with covariates
LD.DRB1.COV <- left_join(LD.DRB1df,LD.COV, by = "LUIN") %>% na.omit() # join together so now we have lt data where each indv matched to their hla allele.

# sort .w.qc col
LD.DRB1.w.qc.df <- as.data.frame(LD.DRB1.w.qc) %>% rownames_to_column(var = "LUIN") # we want the pps in the same format (i.e., matched to lt cloz data)
LD.DRB1.w.qc.COV <- as.data.frame(left_join(LD.COV, LD.DRB1.w.qc.df, by = "LUIN") %>% na.omit() %>% dplyr::select(-c(1:9))) # merge so entries replicated when N has multiple clozuk assays/data points, then get rid of clozuk data - this needs to have same n rows as the cov data so it matches up for the regs. 
LD.DRB1.w.qc.std <- as.data.frame(datawizard::standardise(LD.DRB1.w.qc.COV)) # standardise


# Test linear regression models
#lmer(LD.DRB1.COV$neut_num ~ LD.DRB1.COV$dailydose + LD.DRB1.COV$clozapine + LD.DRB1.COV$norclozapine + LD.DRB1.w.qc.std$ + LD.DRB1.COV$age_at_assay + LD.DRB1.COV$std_age2  + LD.DRB1.COV$SEX + LD.DRB1.COV$TDS + (1|LD.DRB1.COV$LUIN))

DRB1_lmm <- lapply(1:ncol(LD.DRB1.w.qc.std), function(x) lmer(LD.DRB1.COV$neut_num ~ LD.DRB1.COV$dailydose + LD.DRB1.COV$clozapine + LD.DRB1.COV$norclozapine + LD.DRB1.w.qc.std[,x] + LD.DRB1.COV$age_at_assay + LD.DRB1.COV$std_age2  + LD.DRB1.COV$SEX + LD.DRB1.COV$TDS + (1|LD.DRB1.COV$LUIN)))
DRB1_summaries <- lapply(DRB1_lmm, summary)
#DRB1_summaries_b <- lapply(DRB1_summaries, function(x) x$coefficients[2, c(1,4)]) # None survive Bonferroni

# get the nice info out so dont need to faff with the big summaries list 

DRB1_lmm_brief <- NULL
for (i in 1:ncol(LD.DRB1.w.qc)){
  B <- DRB1_summaries[[i]]$coefficients[5,1]
  SE <- DRB1_summaries[[i]]$coefficients[5,2]  
  P <- DRB1_summaries[[i]]$coefficients[5,5]
  allele <- colnames(LD.DRB1.w.qc.std)[i]
  DRB1_lmm_brief <- rbind(DRB1_lmm_brief, data.frame(Gene = "DRB1", HLA_Allele=allele, Estimate=B, Std_Error = SE, P=P))
}



beepr::beep(sound = 2)








#### final bits ####

classic <- c("A","B", "C", "DPB1", "DQA1", "DQB1", "DRB1")

all_excl <- list(excl_A, excl_B, excl_C, excl_DPB1, excl_DQA1, excl_DQB1, excl_DRB1)
for (i in 1:length(classic)){
    print(paste("The following ", length(all_excl[[i]]) ,classic[i] ," alleles were excluded from analysis: "), quote = F)
    print(noquote(all_excl[i]))
    print(' ', quote = F)
}


# make big df
all_lmm_brief <- rbind(A_lmm_brief, B_lmm_brief, C_lmm_brief, DPB1_lmm_brief, DQA1_lmm_brief, DQB1_lmm_brief, DRB1_lmm_brief)

all_lmm_brief$fdr_adj <- p.adjust(p = all_lmm_brief$P, method = "fdr") 
all_lmm_brief$bonf_adj <- p.adjust(p = all_lmm_brief$P, method = "bonferroni") 

all_nom_sig <- NULL
all_fdr_sig <- NULL
all_bonf_sig <- NULL

for (i in 1:length(classic)){
  gene <- classic[i]
  res <- all_lmm_brief %>% filter(Gene == gene)
  nom_sig <- filter(res, P < 0.05)
  if (dim(nom_sig)[1] == 0) { 
    print(paste("No alleles in ", classic[i] ," are significantly associated with ANC"))
  }else{
    print(paste(nrow(nom_sig), " allele(s) in ", classic[i] ," is/are significantly associated with ANC"))
    all_nom_sig <- rbind(all_nom_sig, nom_sig)
  }
  
  fdr_sig <- filter(res, fdr_adj < 0.05)
  if (dim(fdr_sig)[1] == 0) { 
    print(paste("No alleles in ", classic[i], " are significantly associated with ANC (FDR-corrected)"))
  }else {
    print(paste(nrow(fdr_sig), " allele(s) in ", classic[i] ," is/are significantly associated with ANC (FDR-corrected)"))
    all_fdr_sig <- rbind(all_fdr_sig, fdr_sig)
  }
  
  bonf_sig <- filter(res, bonf_adj < 0.05)
  if (dim(fdr_sig)[1] == 0) { 
    print(paste("No alleles in ", classic[i], " are significantly associated with ANC (Bonferroni-corrected)"))
  }else {
    print(paste(nrow(fdr_sig), " allele(s) in ", classic[i] ," is/are significantly associated with ANC (Bonferroni-corrected)"))
    all_fdr_sig <- rbind(all_fdr_sig, fdr_sig)
  }
}


#### end ####

beepr::beep(sound = 3)
