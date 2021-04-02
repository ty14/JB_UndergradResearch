
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
#getting gene information
genes <-genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "ENTREZID"), 
                        keytype='ENSEMBL')
head(genes)

#adding to deg list 
do$genes <- genes
do
do1$genes <- genes
#matrix 

design <- model.matrix(~0+do1$samples$condition4 + do1$samples$batch)
colnames(design) <- gsub("samples", "", colnames(design))
design
colnames(design) <- c("DD", "DS13", "DS14", "batch")  
design
contr.matrix <- makeContrasts(
  DDvsDS14 = DD-DS14, 
  DDvsDS13 = DD-DS13,
  DS13vsDS14 = DS13-DS14,
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(do1, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")


summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

#trouble shooting

design <- model.matrix(~0+do1$samples$condition4)
colnames(design) <- gsub("samples", "", colnames(design))
design
colnames(design) <- c("DD", "DS13", "DS14")  
design
contr.matrix <- makeContrasts(
  DDvsDS14 = DD-DS14, 
  DDvsDS13 = DD-DS13,
  DS13vsDS14 = DS13-DS14,
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(do1, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")


summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))









