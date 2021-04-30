
# libraries 
library(limma)
library(Glimma)
library(edgeR)
library(AnnotationDbi)
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)
library(biomaRt)
library(DESeq2)
library(tidyverse)
library(Mus.musculus)
library(tidyr)
library(gplots)
library(dplyr)
#Getting metadata ready 
coldata <- read_csv("BrainData/raw_data/sample_table.csv")
head(coldata)
str(coldata)

ph <- read_csv("BrainData/raw_data/physiologyData.csv")
head(ph)

ph <- ph %>% 
  unique(.) %>%
  group_by(pre_idbatchcage) %>% 
  pivot_wider(values_from = mean_con_ng_ul, names_from = period) %>% 
  mutate(diff = Post - Pre) %>% ungroup()

ph <- ph %>% filter(time == "1 hr") %>% dplyr::select(1:5,15:17)

coldata <- coldata %>% full_join(ph) %>% unique(.)

#changing things to factors 
coldata$group <- as.factor(coldata$group)
coldata$groupEX <- coldata$group
coldata$Prerank <- as.factor(coldata$Prerank)
coldata$Postrank <- as.factor(coldata$Postrank)
coldata <- coldata %>%  dplyr::select(-group)

##getting condition2
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "HeldRank", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 3, "HeldRank", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 4, "HeldRank", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders"| coldata$condition == "ascenders", "ChangedRank", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control", "Control", coldata$condition1)

#another way to look at condition
coldata$condition2 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "HeldRank", coldata$condition)
coldata$condition2 <- ifelse(coldata$condition == "same" & coldata$Prerank == 3, "HeldRank", coldata$condition2)
coldata$condition2 <- ifelse(coldata$condition == "same" & coldata$Prerank == 4, "HeldRank", coldata$condition2)
coldata$condition2 <- ifelse(coldata$condition == "descenders", "Descenders", coldata$condition2)
coldata$condition2 <- ifelse(coldata$condition == "ascenders", "Ascenders", coldata$condition2)
coldata$condition2 <- ifelse(coldata$condition == "control", "Control", coldata$condition2)

## condition3 because we are more interested in the comparison between 
# Dom-Dom to Descenders (DOM to SUB)
# Sub-Sub to Ascenders (Sub to DOM)
coldata$condition3 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DD", coldata$condition)
coldata$condition3 <- ifelse(coldata$condition == "same" & coldata$Prerank == 3, "SS", coldata$condition3)
coldata$condition3 <- ifelse(coldata$condition == "same" & coldata$Prerank == 4, "SS", coldata$condition3)
coldata$condition3 <- ifelse(coldata$condition == "descenders", "DS", coldata$condition3)
coldata$condition3 <- ifelse(coldata$condition == "ascenders", "SD", coldata$condition3)
coldata$condition3 <- ifelse(coldata$condition == "control", "Control", coldata$condition3)

## condition4 because we are more interested in the comparison between 
# Dom-Dom to Descenders (DOM to SUB) (3->1) (4->1)
# Sub-Sub to Ascenders (Sub to DOM) (1->3) (1->4)
coldata$condition4 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DD", coldata$condition)
coldata$condition4 <- ifelse(coldata$condition == "same" & coldata$Prerank == 3, "SS", coldata$condition4)
coldata$condition4 <- ifelse(coldata$condition == "same" & coldata$Prerank == 4, "SS", coldata$condition4)
coldata$condition4 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DS14", coldata$condition4)
coldata$condition4 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 3 & coldata$Prerank == 1, "DS13", coldata$condition4)
coldata$condition4 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "SD41", coldata$condition4)
coldata$condition4 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 3 & coldata$Postrank == 1, "SD31", coldata$condition4)
coldata$condition4 <- ifelse(coldata$condition == "control", "Control", coldata$condition4)

# coldata$condition4 <- factor(coldata$condition4, level = c("Control", "DD", "SS", "SD41", "SD31", "DS14", "DS13"))

# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

#Pre
coldata <- transform(coldata, Pre_norm = (Pre - min(Pre)) / (max(Pre) - min(Pre)))
#Post
coldata <- transform(coldata, Post_norm = (Post - min(Post)) / (max(Post) - min(Post)))
#Difference in cort 
coldata <- transform(coldata, diff_norm = (diff - min(diff)) / (max(diff) - min(diff)))

# First PFC data 
pfc_data<- coldata %>% filter(region == "PFC")%>% arrange(condition1) %>% arrange(groupEX) 

row.names <- pfc_data$SampleName

pfc_data1 <- pfc_data%>% dplyr::select(-SampleName)

row.names(pfc_data1) <- row.names #Assigning row names from as sample names  
head(pfc_data1)

#bring in physiology data 
ph <- read_csv("BrainData/raw_data/physiologyData.csv")
head(ph)

# Bring in count data for mPFC
p_countdata <- read_csv("BrainData/raw_data/PFC_counts.csv")
p_countdata$X -> nrows
p_count <- p_countdata[,-1]
row.names(p_count) <- nrows
p_count <-as.data.frame(p_count)

#checks
all (row.names(pfc_data1) %in% colnames(p_count)) #check 

p_count <- p_count[,rownames(pfc_data1)]
all(rownames(pfc_data1) == colnames(p_count)) #check

# Create DGEList object
do <- DGEList(p_count)
do


#getting sample data in list
do$samples

samplenames<- substring(colnames(do),1, nchar(colnames(do)))
samplenames
colnames(do) <- samplenames


do$samples <-cbind(SampleName=rownames(do$samples), do$samples)

do$samples <- do$samples %>% full_join(pfc_data)
rownames(do$samples) <- do$samples$SampleName

do$samples <- do$samples %>% dplyr::select(-SampleName)
#That was a pain in the butt

#filter out controls 1. just look at  DD DS  2. just look at SS and SD (1->3) (1->4) and vice versa 

do<- do[,do$samples$condition4 %in% c("DD","SS", "DS14", "DS13", "SD41", "SD31")]
do$samples$condition4


do1<- do[,do$samples$condition4 %in% c("DD","DS14", "DS13")] 
do1$samples$condition4

do3<- do[,do$samples$condition4 %in% c("SS","SD41", "SD31")]
do3$samples$condition4





geneid <- rownames(do$counts)
geneid <- rownames(do1$counts)
geneid <- rownames(do3$counts)
#getting gene information
genes <-genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "ENTREZID"), 
                        keytype='ENSEMBL')
head(genes)

#adding to deg list 
do$genes <- genes
do
do1$genes <- genes
do3$genes <- genes

# remove FALSE rows 
do1 <- do1[ which(rownames(do1$counts) %in% do1$genes$ENSEMBL),]
rownames(do1$counts) %in% do1$genes$ENSEMBL

#data processing do1
mean(do1$samples$lib.size)
cpm <- cpm(do1)
lcpm <- cpm(do1, log=TRUE)

L <- mean(do1$samples$lib.size) * 1e-6
M <- median(do1$samples$lib.size) * 1e-6
c(L, M)

summary(lcpm)

table(rowSums(do1$counts==0)==21)
keep.exprs <- filterByExpr(do1, samples= do1$samples$condition4)
do1 <- do1[keep.exprs,, keep.lib.sizes=FALSE]
dim(do1)

lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(do1)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(do1, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

do1 <- calcNormFactors(do1, method = "TMM")
do1$samples$norm.factors

do1 <- do1
do1$samples$norm.factors <- 1
do1$counts[,1] <- ceiling(do1$counts[,1]*0.05)
do1$counts[,2] <- do1$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(do1, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
do1 <- calcNormFactors(do1)  
do1$samples$norm.factors

lcpm <- cpm(do1, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

#matrix do1

design <- model.matrix(~0+do1$samples$condition4)
colnames(design) <- gsub("samples", "", colnames(design))
design
colnames(design) <- c("DD", "DS13", "DS14")  
design
contr.matrix <- makeContrasts(
  DS14vsDD = DS14-DD, 
  DS13vsDD = DS13-DD,
  DS14vsDS13 = DS14-DS13,
  levels = colnames(design))
contr.matrix

# How many random sampling
R = 5000

set.seed(312)

p.dl.rand = list()

for(i in 1 : R){
  print(paste("Starting on Permutation", i))
  
  # Randomize the traits
  counts.rand = sample(do1$counts)
  # Model
  design.dl.rand = model.matrix(~0 + do1$samples$condition4 + counts.rand)
  colnames(design.dl.rand) = c("DS14vsDD", "DS13vsDD", "DS14vsDS13", "count")
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(do1, design.dl.rand, plot=F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  vfit.dl.rand = contrasts.fit(vfit.dl.rand, contrasts = contr.matrix)
  efit.dl.rand = eBayes(vfit.dl.rand)
  p.dl.rand[[i]] = efit.dl.rand[["p.value"]]
}

q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

for(i in 1 : R){
  print(paste("Calculating Permutation", i))
  
  temp = p.dl.rand[[i]]
  
  for(c in 1 : 4){
    for(r in 1 : nrow(p.dl.limma)){
      if(temp[r, c] <= p.dl.limma[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
colnames(q.dl) = c("DS41vsDD", "DS31vsDD", "DS41vsDS31", "count")
q.dl = as.data.frame(q.dl)

#Venn Diagram 

par(mfrow=c(1,2))
v <- voom(do1, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

dtp <- decideTests(efit, adjust.method="none", p.value=0.05,lfc=0.2)
summary(dtp)
de.common <- which(dtp[,1]!=0 & dtp[,2]!=0)
length(de.common)
head(efit$genes$SYMBOL[de.common], n=20)
vennDiagram(dtp[,1:2], circle.col=c("turquoise", "salmon"))
write.fit(efit, dtp, file="results.txt")
# Examining Individual DE
DS14.vs.DD <- topTable(efit, coef=1, n=Inf)
DS14DD <- DS14.vs.DD %>% rownames_to_column(., var = "ENSEMBLS") %>% as_tibble(.) %>% mutate(contrast = "DS14 - DD") %>% full_join(genes)
DS13.vs.DD <- topTable(efit, coef=2, n=Inf)
DS13DD <- DS13.vs.DD %>% rownames_to_column(., var = "ENSEMBLS") %>% as_tibble(.) %>% mutate(contrast = "DS13 - DD") %>% full_join(genes)
DS14.vs.DS13 <- topTable(efit, coef=3, n=Inf)
DS14DS13 <- DS14.vs.DS13 %>% rownames_to_column(., var = "ENSEMBLS") %>% as_tibble(.) %>% mutate(contrast = "DS14 - DS13") %>% full_join(genes)
head(DS14.vs.DD)
head(DS13.vs.DD)
head(DS14.vs.DS13)


#rbind and save
all_genes <- DS14DD %>% rbind(DS13DD,DS14DS13)
head(all_genes)

write.csv(all_genes, "datacleann/Condition4Cort_PFCgenes.csv", row.names = F)


plotMD(efit, column=1, status=dtp[,1], main=colnames(efit)[1], 
       xlim=c(-5,13))

library(gplots)
DS14.vs.DD.topgenes <- DS14.vs.DD$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% DS14.vs.DD.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm, scale="row",
          labRow=v$genes$SYMBOL[i], labCol=do1$samples$condition4, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

# gene set testing 

load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
idx <- ids2indices(Mm.c2,id=rownames(do1$genes))
cam.DS14vsDD <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.DS14vsDD,5)

cam.DS13vsDD <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.DS13vsDD,5)

cam.DS14vsDS13 <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.DS14vsDS13,5)

barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="DS14vsDS13")



#data processing do3
mean(do3$samples$lib.size)
cpm <- cpm(do3)
lcpm <- cpm(do3, log=TRUE)

L <- mean(do3$samples$lib.size) * 1e-6
M <- median(do3$samples$lib.size) * 1e-6
c(L, M)

summary(lcpm)

table(rowSums(do3$counts==0)==27)
keep.exprs <- filterByExpr(do3, samples= do3$samples$condition4)
do3 <- do3[keep.exprs,, keep.lib.sizes=FALSE]
dim(do3)

lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(do3)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(do3, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

do3 <- calcNormFactors(do3, method = "TMM")
do3$samples$norm.factors

do3 <- do3
do3$samples$norm.factors <- 1
do3$counts[,1] <- ceiling(do3$counts[,1]*0.05)
do3$counts[,2] <- do3$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(do3, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
do3 <- calcNormFactors(do3)  
do3$samples$norm.factors

lcpm <- cpm(do3, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

#matrix do3

design <- model.matrix(~0+do3$samples$condition4)
colnames(design) <- gsub("samples", "", colnames(design))
design
colnames(design) <- c("SS", "SD41", "SD31")  
design
contr.matrix <- makeContrasts(
  SD41vsSS = SD41-SS, 
  SD31vsSS = SD31-SS,
  SD41vsSD31 = SD41-SD31,
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(do3, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

dtp <- decideTests(efit, adjust.method="none", p.value=0.05,lfc=0.2)
summary(dtp)
de.common <- which(dtp[,1]!=0 & dtp[,2]!=0)
length(de.common)
head(efit$genes$SYMBOL[de.common], n=20)
vennDiagram(dtp[,1:2], circle.col=c("turquoise", "salmon"))

# Examining Individual DE
SD41.vs.SS <- topTable(efit, coef=1, n=Inf)
SD41SS <- SD41.vs.SS %>% rownames_to_column(., var = "ENSEMBLS") %>% as_tibble(.) %>% mutate(contrast = "SD41 - SS") %>% full_join(genes)
SD31.vs.SS <- topTable(efit, coef=2, n=Inf)
SD31SS <- SD31.vs.SS %>% rownames_to_column(., var = "ENSEMBLS") %>% as_tibble(.) %>% mutate(contrast = "SD31 - SS") %>% full_join(genes)
SD41.vs.SD31 <- topTable(efit, coef=3, n=Inf)
SD41SD31 <- SD41.vs.SD31 %>% rownames_to_column(., var = "ENSEMBLS") %>% as_tibble(.) %>% mutate(contrast = "SD41 - SD31") %>% full_join(genes)
head(SD41.vs.SS)
head(SD31.vs.SS)
head(SD41.vs.SD31)

#rbind and save
all_genes2 <- SD41SS %>% rbind(SD31SS,SD41SD31)
head(all_genes2)

write.csv(all_genes2, "manuscript/brain/LimmaAnalysis/Condition3Cort_PFCgenes.csv", row.names = F)

plotMD(efit, column=1, status=dtp[,1], main=colnames(efit)[1], 
       xlim=c(-5,13))
library(gplots)
SD41.vs.SS.topgenes <- SD41.vs.SS$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% SD41.vs.SS.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm, scale="row",
          labRow=v$genes$SYMBOL[i], labCol=do3$samples$condition4, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
# gene set testing 

load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
idx <- ids2indices(Mm.c2,id=rownames(v$genes))
cam.SD31vsSS <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.SD31vsSS,5)

cam.SD41vsSS <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.SD41vsSS,5)

cam.SD41vsSD31 <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.SD41vsSD31,5)

barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="SD41vsSD31")
