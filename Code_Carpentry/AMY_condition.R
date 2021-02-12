# LIMMA ANALYSIS

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

# coldata$condition3 <- factor(coldata$condition3, level = c("Control", "DD", "SS", "SD", "DS"))

# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

#Pre
coldata <- transform(coldata, Pre_norm = (Pre - min(Pre)) / (max(Pre) - min(Pre)))
#Post
coldata <- transform(coldata, Post_norm = (Post - min(Post)) / (max(Post) - min(Post)))
#Difference in cort 
coldata <- transform(coldata, diff_norm = (diff - min(diff)) / (max(diff) - min(diff)))

# First AMY data 
amy_data<- coldata %>% filter(region == "AMY")%>% arrange(condition1) %>% arrange(groupEX) 

row.names <- amy_data$SampleName

amy_data1 <- amy_data%>% dplyr::select(-SampleName)

row.names(amy_data1) <- row.names #Assigning row names from as sample names  
head(amy_data1)

#bring in physiology data 
ph <- read_csv("manuscript/cort/physiologyData.csv")
head(ph)

# Bring in count data for mAMY
a_countdata <- read_csv("manuscript/brain/AMY_counts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
a_count <-as.data.frame(a_count)

#checks
all (row.names(amy_data1) %in% colnames(a_count)) #check 

a_count <- a_count[,rownames(amy_data1)]
all(rownames(amy_data1) == colnames(a_count)) #check

# Create DGEList object
do <- DGEList(a_count)
do


#getting sample data in list
do$samples

samplenames<- substring(colnames(do),1, nchar(colnames(do)))
samplenames
colnames(do) <- samplenames


do$samples <-cbind(SampleName=rownames(do$samples), do$samples)

do$samples <- do$samples %>% full_join(amy_data)
rownames(do$samples) <- do$samples$SampleName

do$samples <- do$samples %>% dplyr::select(-SampleName)
#That was a pain in the butt

#filter out controls 1. just look at  DD DS  2. just look at SS and SD

do<- do[,do$samples$condition3 %in% c("DD","SS", "DS", "SD")]
do$samples$condition3


do1<- do[,do$samples$condition3 %in% c("DD","DS")]
do1$samples$condition3

do2<- do[,do$samples$condition3 %in% c("SS","SD")]
do2$samples$condition3



geneid <- rownames(do$counts)
#getting gene information
genes <-genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "ENTREZID"), 
                        keytype='ENSEMBL')
head(genes)

#adding to deg list 
do$genes <- genes
do



